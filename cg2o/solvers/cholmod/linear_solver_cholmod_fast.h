/*
Contributed to cg2o; distributed under the same BSD 3-Clause license as the
rest of the library (see the header of linear_solver_eigen_umfpack_lu.h).
*/

#ifndef CG2O_LINEAR_SOLVER_CHOLMOD_FAST_H
#define CG2O_LINEAR_SOLVER_CHOLMOD_FAST_H

#include <cholmod.h>

#include <algorithm>
#include <cmath>
#include <functional>
#include <numeric>
#include <utility>
#include <vector>

#include "g2o/core/batch_stats.h"
#include "g2o/core/linear_solver.h"
#include "g2o/core/marginal_covariance_cholesky.h"
#include "g2o/core/matrix_structure.h"
#include "g2o/stuff/logger.h"
#include "g2o/stuff/timeutil.h"

namespace cg2o {

/**
 * \brief CHOLMOD (SuiteSparse) for any g2o system, without the per-call
 * overhead of g2o::LinearSolverCholmod.
 *
 * Differences to g2o::LinearSolverCholmod:
 *  - the fill-reducing ordering and the symbolic analysis are kept across
 *    optimize() calls while the sparsity pattern is unchanged (g2o frees the
 *    factor in init(), i.e. at every optimize(); in MPC the pattern never
 *    changes);
 *  - pre-permutation (Options::prepermute, default on): the ordering P is
 *    computed once, and the values g2o assembles are scattered directly into
 *    the upper triangle of P A P^T through a precomputed map. CHOLMOD then
 *    factors with the natural ordering and skips the permuted, transposed
 *    copy of the matrix it otherwise builds inside every cholmod_factorize;
 *  - factor type (Options::supernodal): CHOLMOD_AUTO (default) lets CHOLMOD
 *    pick simplicial LDL^T for sparse factors and supernodal LL^T for denser
 *    ones (faster on large SPD problems). Supernodal LL^T cannot factor an
 *    indefinite matrix; if it fails, the solver switches to simplicial LDL^T
 *    for good (one warning). Simplicial LDL^T factors SPD matrices and also
 *    symmetric quasi-definite KKT matrices ([H J^T; J -sigma I], e.g. cg2o
 *    ISPD/BIPM with sigma in the equality edges, or Options::delta > 0).
 *    CHOLMOD_SIMPLICIAL skips the automatic choice (KKT systems);
 *    CHOLMOD_SUPERNODAL forces LL^T (SPD only, no fallback);
 *  - cholmod_solve2 reuses its work vectors instead of allocating a new
 *    result vector on every solve;
 *  - optional quasi-definite mode (Options::delta > 0): -delta on exactly-zero
 *    diagonal entries with adaptive refinement against the unregularized
 *    matrix (only when some entry was regularized), as in
 *    LinearSolverEigenQDLDLT;
 *  - marginal covariances for SPD systems.
 * g2o's block ordering (setBlockOrdering(true), the LinearSolverCCS default)
 * computes AMD on the block structure, like g2o::LinearSolverCholmod;
 * setBlockOrdering(false) uses Options::ordering on the scalar matrix.
 */
template <typename MatrixType>
class LinearSolverCholmodFast : public g2o::LinearSolverCCS<MatrixType> {
public:
  struct Options {
    int supernodal = CHOLMOD_AUTO; // or CHOLMOD_SIMPLICIAL (LDL^T),
                                   // CHOLMOD_SUPERNODAL (LL^T, SPD only)
    bool prepermute =
        true; // scatter into P A P^T (no CHOLMOD-internal permuted copy)
    int ordering = CHOLMOD_AMD; // scalar ordering if blockOrdering() is false
                                // (AMD, METIS, NESDIS, ...)
    double delta = 0.0;         // > 0 enables quasi-definite mode
    int refinementSteps =
        5; // max own refinement steps (only if entries were regularized)
    double refinementTol =
        1e-12; // stop when max|r| <= tol * max|b| (or stalled)
  };

  LinearSolverCholmodFast() : g2o::LinearSolverCCS<MatrixType>() {
    cholmod_start(&_c);
  }
  ~LinearSolverCholmodFast() override {
    release();
    cholmod_free_dense(&_X, &_c);
    cholmod_free_dense(&_Y, &_c);
    cholmod_free_dense(&_E, &_c);
    cholmod_finish(&_c);
  }
  LinearSolverCholmodFast(const LinearSolverCholmodFast &) = delete;
  LinearSolverCholmodFast &operator=(const LinearSolverCholmodFast &) = delete;

  //! change before the first solve (a change of supernodal/prepermute/ordering
  //! triggers a new analysis at the next optimize())
  Options &options() { return _opt; }
  //! diagonal entries that received -delta in the last factorization
  int lastRegularizedCount() const { return _numRegularized; }
  //! refinement steps taken in the last solve
  int lastRefinementSteps() const { return _lastRefinementSteps; }
  //! true if CHOLMOD chose (or was told to use) the supernodal LL^T factor
  bool isSupernodal() const { return _L != nullptr && _L->is_super; }

  bool init() override {
    _init = true; // pattern is re-checked; analysis kept if unchanged
    return true;
  }

  // g2o's GN/LM apply the update vector even when solve() fails, so a
  // failed solve leaves x = 0 (no step) instead of stale values.
  bool solve(const g2o::SparseBlockMatrix<MatrixType> &A, double *x,
             double *b) override {
    if (solveImpl(A, x, b))
      return true;
    std::fill(x, x + A.rows(), 0.0);
    return false;
  }

protected:
  bool solveImpl(const g2o::SparseBlockMatrix<MatrixType> &A, double *x,
                 double *b) {
    if (!factorize(A))
      return false;
    const double t = g2o::get_monotonic_time();
    for (int k = 0; k < _n; ++k)
      _bw[k] = b[_perm[k]]; // into factor order
    if (!solveFactored(_bw.data(), _xw.data()))
      return false;
    _lastRefinementSteps = 0;
    if (_numRegularized > 0)
      refine();
    for (int k = 0; k < _n; ++k)
      x[_perm[k]] = _xw[k]; // back to g2o order
    for (int k = 0; k < _n; ++k)
      if (!std::isfinite(x[k])) {
        G2O_ERROR("LinearSolverCholmodFast: solve produced non-finite values");
        return false;
      }
    g2o::G2OBatchStatistics *s = g2o::G2OBatchStatistics::globalStats();
    if (s) {
      s->timeNumericDecomposition += g2o::get_monotonic_time() - t;
      s->choleskyNNZ = static_cast<size_t>(_lnz);
    }
    return true;
  }

  Options _opt;
  cholmod_common _c;
  cholmod_factor *_L = nullptr;
  cholmod_dense *_X = nullptr, *_Y = nullptr, *_E = nullptr; // solve2 work
  bool _init = true;
  bool _forceSimplicial = false; // set when supernodal LL^T failed (AUTO)
  int _n = 0;
  double _lnz = 0;
  // settings the current analysis was made with
  int _analyzedMode = -1, _analyzedOrdering = -1;
  bool _analyzedPrepermute = false, _analyzedBlockOrdering = false;
  // g2o's upper CCS (source) and the matrix CHOLMOD factors (upper CCS of
  // P A P^T with prepermute, else the same pattern as the source)
  std::vector<int> _Cp, _Ci, _Wp, _Wi, _map, _diagPos, _perm;
  std::vector<double> _Cx, _Wx, _bw, _xw, _r, _dx;
  std::vector<char> _regMask;
  int _numRegularized = 0;
  int _lastRefinementSteps = 0;
  g2o::MatrixStructure _blockStructure;

  void release() {
    if (_L)
      cholmod_free_factor(&_L, &_c);
  }

  static cholmod_sparse view(int n, std::vector<int> &p, std::vector<int> &i,
                             double *x, int xtype) {
    cholmod_sparse s;
    s.nrow = s.ncol = n;
    s.nzmax = i.size();
    s.p = p.data();
    s.i = i.data();
    s.nz = nullptr;
    s.x = x;
    s.z = nullptr;
    s.stype = 1; // upper triangle of a symmetric matrix
    s.itype = CHOLMOD_INT;
    s.xtype = xtype;
    s.dtype = CHOLMOD_DOUBLE;
    s.sorted = 1;
    s.packed = 1;
    return s;
  }

  void setStat(bool symbolic, double dt) {
    g2o::G2OBatchStatistics *s = g2o::G2OBatchStatistics::globalStats();
    if (!s)
      return;
    if (symbolic)
      s->timeSymbolicDecomposition = dt;
    else
      s->timeNumericDecomposition = dt;
  }

  int currentMode() const {
    return _forceSimplicial ? CHOLMOD_SIMPLICIAL : _opt.supernodal;
  }

  bool settingsChanged() const {
    return _analyzedMode != currentMode() ||
           _analyzedPrepermute != _opt.prepermute ||
           _analyzedOrdering != _opt.ordering ||
           _analyzedBlockOrdering != this->blockOrdering();
  }

  bool factorize(const g2o::SparseBlockMatrix<MatrixType> &A) {
    if (_init) {
      const double t0 = g2o::get_monotonic_time();
      this->initMatrixStructure(A);
      const int n = A.rows(), nz = A.nonZeros();
      std::vector<int> Cp(n + 1), Ci(nz);
      _Cx.resize(nz);
      // A.nonZeros() counts whole blocks; the upper triangle has fewer entries
      const int nzU =
          this->_ccsMatrix->fillCCS(Cp.data(), Ci.data(), _Cx.data(), true);
      Ci.resize(nzU);
      _Cx.resize(nzU);
      const bool same = _L != nullptr && n == _n && Cp == _Cp && Ci == _Ci &&
                        !settingsChanged();
      _init = false;
      if (!same) {
        release();
        _n = n;
        _Cp.swap(Cp);
        _Ci.swap(Ci);
        if (!analyze(A))
          return false;
        setStat(true, g2o::get_monotonic_time() - t0);
      }
    } else {
      this->_ccsMatrix->fillCCS(_Cx.data(), true);
    }
    const double t1 = g2o::get_monotonic_time();
    if (!factorizeValues()) {
      // CHOLMOD_AUTO picked supernodal LL^T, which cannot factor an
      // indefinite (KKT) matrix: switch to simplicial LDL^T for good
      if (_L != nullptr && _L->is_super &&
          _opt.supernodal != CHOLMOD_SUPERNODAL && !_forceSimplicial) {
        G2O_WARN("LinearSolverCholmodFast: supernodal LL^T failed (indefinite "
                 "matrix?), switching to simplicial LDL^T");
        _forceSimplicial = true;
        release();
        if (!analyze(A) || !factorizeValues())
          return failed(A);
      } else {
        return failed(A);
      }
    }
    setStat(false, g2o::get_monotonic_time() - t1);
    return true;
  }

  bool failed(const g2o::SparseBlockMatrix<MatrixType> &A) {
    if (this->writeDebug()) {
      G2O_ERROR(
          "LinearSolverCholmodFast: factorization failed, writing debug.txt");
      A.writeOctave("debug.txt");
    } else {
      G2O_DEBUG("LinearSolverCholmodFast: factorization failed");
    }
    return false;
  }

  // scatter the values into the factored matrix, regularize, factorize
  bool factorizeValues() {
    const int nzU = static_cast<int>(_Cx.size());
    for (int e = 0; e < nzU; ++e)
      _Wx[_map[e]] = _Cx[e];
    _numRegularized = 0;
    if (_opt.delta > 0.0)
      for (int c = 0; c < _n; ++c) {
        _regMask[c] = _diagPos[c] >= 0 && _Wx[_diagPos[c]] == 0.0;
        if (_regMask[c]) {
          _Wx[_diagPos[c]] = -_opt.delta;
          ++_numRegularized;
        }
      }
    cholmod_sparse W = view(_n, _Wp, _Wi, _Wx.data(), CHOLMOD_REAL);
    cholmod_factorize(&W, _L, &_c);
    // minor < n: a pivot was zero (LDL^T) or not positive (LL^T)
    return _c.status >= CHOLMOD_OK && _L->minor >= _L->n;
  }

  // ordering + symbolic analysis for the current pattern and settings
  bool analyze(const g2o::SparseBlockMatrix<MatrixType> &A) {
    const int n = _n, nzU = static_cast<int>(_Ci.size());
    const int mode = currentMode();
    cholmod_sparse S = view(n, _Cp, _Ci, nullptr, CHOLMOD_PATTERN);
    std::vector<int> given; // block ordering expanded to the scalar matrix
    if (this->blockOrdering()) {
      A.fillBlockStructure(_blockStructure);
      cholmod_sparse B;
      B.nrow = B.ncol = _blockStructure.n;
      B.nzmax = _blockStructure.nzMax();
      B.p = _blockStructure.Ap;
      B.i = _blockStructure.Aii;
      B.nz = nullptr;
      B.x = B.z = nullptr;
      B.stype = 1;
      B.itype = CHOLMOD_INT;
      B.xtype = CHOLMOD_PATTERN;
      B.dtype = CHOLMOD_DOUBLE;
      B.sorted = 1;
      B.packed = 1;
      Eigen::VectorXi blockPerm(_blockStructure.n), scalarPerm(n);
      if (!cholmod_amd(&B, nullptr, 0, blockPerm.data(), &_c)) {
        G2O_ERROR("LinearSolverCholmodFast: block AMD failed");
        return false;
      }
      this->blockToScalarPermutation(A, blockPerm, scalarPerm);
      given.assign(scalarPerm.data(), scalarPerm.data() + n);
    }
    auto runAnalysis = [&](cholmod_sparse *M, int *perm, int ordering,
                           int postorder, int supernodal) {
      _c.nmethods = 1;
      _c.method[0].ordering = perm ? CHOLMOD_GIVEN : ordering;
      _c.postorder = postorder;
      _c.supernodal = supernodal;
      return perm ? cholmod_analyze_p(M, perm, nullptr, 0, &_c)
                  : cholmod_analyze(M, &_c);
    };
    _perm.resize(n);
    if (_opt.prepermute) {
      // 1. fill-reducing ordering (incl. postorder) of the source pattern
      cholmod_factor *L0 =
          runAnalysis(&S, given.empty() ? nullptr : given.data(), _opt.ordering,
                      1, CHOLMOD_SIMPLICIAL);
      if (!L0) {
        G2O_ERROR("LinearSolverCholmodFast: ordering failed");
        return false;
      }
      std::copy((int *)L0->Perm, (int *)L0->Perm + n, _perm.begin());
      cholmod_free_factor(&L0, &_c);
      // 2. pattern of the upper triangle of P A P^T and the scatter map
      std::vector<int> pinv(n);
      for (int k = 0; k < n; ++k)
        pinv[_perm[k]] = k;
      std::vector<std::vector<std::pair<int, int>>> col(
          n); // (row, source entry)
      for (int j = 0; j < n; ++j)
        for (int k = _Cp[j]; k < _Cp[j + 1]; ++k) {
          const int a = pinv[_Ci[k]], b = pinv[j];
          col[std::max(a, b)].push_back({std::min(a, b), k});
        }
      _Wp.assign(n + 1, 0);
      _Wi.resize(nzU);
      _map.assign(nzU, -1);
      _diagPos.assign(n, -1);
      int p = 0;
      for (int j = 0; j < n; ++j) {
        std::sort(col[j].begin(), col[j].end());
        _Wp[j] = p;
        for (auto &re : col[j]) {
          _Wi[p] = re.first;
          _map[re.second] = p;
          if (re.first == j)
            _diagPos[j] = p;
          ++p;
        }
      }
      _Wp[n] = p;
      // 3. symbolic analysis of P A P^T with the natural ordering
      cholmod_sparse W = view(n, _Wp, _Wi, nullptr, CHOLMOD_PATTERN);
      _L = runAnalysis(&W, nullptr, CHOLMOD_NATURAL, 0, mode);
    } else {
      std::iota(_perm.begin(), _perm.end(), 0);
      _Wp = _Cp;
      _Wi = _Ci;
      _map.resize(nzU);
      std::iota(_map.begin(), _map.end(), 0);
      _diagPos.assign(n, -1);
      for (int j = 0; j < n; ++j)
        for (int k = _Cp[j]; k < _Cp[j + 1]; ++k)
          if (_Ci[k] == j)
            _diagPos[j] = k;
      _L = runAnalysis(&S, given.empty() ? nullptr : given.data(),
                       _opt.ordering, 1, mode);
    }
    if (!_L) {
      G2O_ERROR("LinearSolverCholmodFast: symbolic analysis failed ({})",
                _c.status);
      return false;
    }
    _lnz = _c.lnz;
    _Wx.assign(nzU, 0.0);
    _bw.resize(n);
    _xw.resize(n);
    _r.resize(n);
    _dx.resize(n);
    _regMask.assign(n, 0);
    _analyzedMode = mode;
    _analyzedPrepermute = _opt.prepermute;
    _analyzedOrdering = _opt.ordering;
    _analyzedBlockOrdering = this->blockOrdering();
    return true;
  }

  // solve with the current factor; b and x in factor (W) order
  bool solveFactored(const double *b, double *x) {
    cholmod_dense B;
    B.nrow = B.d = _n;
    B.ncol = 1;
    B.nzmax = _n;
    B.x = const_cast<double *>(b);
    B.z = nullptr;
    B.xtype = CHOLMOD_REAL;
    B.dtype = CHOLMOD_DOUBLE;
    if (!cholmod_solve2(CHOLMOD_A, _L, &B, nullptr, &_X, nullptr, &_Y, &_E,
                        &_c)) {
      G2O_ERROR("LinearSolverCholmodFast: cholmod_solve2 failed ({})",
                _c.status);
      return false;
    }
    std::copy((double *)_X->x, (double *)_X->x + _n, x);
    return true;
  }

  // adaptive refinement against the unregularized matrix (see QDLDLT)
  void refine() {
    double normB = 0.0, prev = 1e300;
    for (int i = 0; i < _n; ++i)
      normB = std::max(normB, std::abs(_bw[i]));
    for (int it = 0; it < _opt.refinementSteps; ++it) {
      std::copy(_bw.begin(), _bw.end(), _r.begin());
      for (int j = 0; j < _n; ++j) // symmetric product from the upper CCS
        for (int k = _Wp[j]; k < _Wp[j + 1]; ++k) {
          const int i = _Wi[k];
          const double v = (k == _diagPos[j] && _regMask[j]) ? 0.0 : _Wx[k];
          _r[i] -= v * _xw[j];
          if (i != j)
            _r[j] -= v * _xw[i];
        }
      double rn = 0.0;
      for (int i = 0; i < _n; ++i)
        rn = std::max(rn, std::abs(_r[i]));
      if (rn <= _opt.refinementTol * normB || rn > 0.5 * prev)
        break;
      prev = rn;
      if (!solveFactored(_r.data(), _dx.data()))
        break;
      for (int i = 0; i < _n; ++i)
        _xw[i] += _dx[i];
      ++_lastRefinementSteps;
    }
  }

  /**
   * Marginal covariances for SPD systems: a copy of the factor is converted
   * to the simplicial LL^T form g2o's MarginalCovarianceCholesky expects.
   */
  bool solveBlocks_impl(
      const g2o::SparseBlockMatrix<MatrixType> &A,
      std::function<void(g2o::MarginalCovarianceCholesky &)> compute) override {
    if (!factorize(A))
      return false;
    if (_numRegularized > 0) {
      G2O_WARN("LinearSolverCholmodFast: marginal covariances need a positive "
               "definite system");
      return false;
    }
    if (!_L->is_ll && !_L->is_super) { // simplicial LDL^T: D on the diagonal
      const int *Lp = static_cast<const int *>(_L->p);
      const double *Lx = static_cast<const double *>(_L->x);
      for (int j = 0; j < _n; ++j)
        if (!(Lx[Lp[j]] > 0.0)) {
          G2O_WARN("LinearSolverCholmodFast: marginal covariances need a "
                   "positive definite system");
          return false;
        }
    }
    cholmod_factor *F = cholmod_copy_factor(_L, &_c);
    if (!F)
      return false;
    // simplicial, LL^T, packed, monotonic; first entry of each column =
    // diagonal
    if (!cholmod_change_factor(CHOLMOD_REAL, 1, 0, 1, 1, F, &_c)) {
      cholmod_free_factor(&F, &_c);
      return false;
    }
    // overall ordering: factor row k is original variable _perm[F->Perm[k]]
    const int *fp = static_cast<const int *>(F->Perm);
    std::vector<int> pinv(_n);
    for (int k = 0; k < _n; ++k)
      pinv[_perm[fp ? fp[k] : k]] = k;
    g2o::MarginalCovarianceCholesky mcc;
    mcc.setCholeskyFactor(_n, static_cast<int *>(F->p),
                          static_cast<int *>(F->i), static_cast<double *>(F->x),
                          pinv.data());
    compute(mcc);
    cholmod_free_factor(&F, &_c);
    g2o::G2OBatchStatistics *s = g2o::G2OBatchStatistics::globalStats();
    if (s)
      s->choleskyNNZ = static_cast<size_t>(_lnz);
    return true;
  }
};

} // namespace cg2o

#endif // CG2O_LINEAR_SOLVER_CHOLMOD_FAST_H
