/*
Contributed to cg2o; distributed under the same BSD 3-Clause license as the
rest of the library (see the header of linear_solver_eigen_umfpack_lu.h).
*/

#ifndef CG2O_LINEAR_SOLVER_UMFPACK_FAST_H
#define CG2O_LINEAR_SOLVER_UMFPACK_FAST_H

#include <umfpack.h>

#include <algorithm>
#include <cmath>
#include <functional>
#include <vector>

#include "g2o/core/batch_stats.h"
#include "g2o/core/linear_solver.h"
#include "g2o/core/marginal_covariance_cholesky.h"
#include "g2o/stuff/logger.h"
#include "g2o/stuff/timeutil.h"

namespace cg2o {

/**
 * \brief UMFPACK (sparse LU with pivoting) for any g2o system, tuned for the
 * symmetric matrices g2o assembles (SPD or symmetric indefinite KKT).
 *
 * Differences to LinearSolverEigenUmfPackLU:
 *  - calls the UMFPACK C API directly; the full (both triangles) CSC pattern is
 *    built once and new values are scattered into it, instead of building a
 *    full symmetric copy through Eigen on every factorization;
 *  - symbolic analysis is kept across optimize() calls while the sparsity
 *    pattern is unchanged;
 *  - the symmetric strategy is used (AMD on A + A^T, diagonal pivots
 *    preferred), and UMFPACK's built-in iterative refinement is configurable
 *    (default off: it costs extra solves and products every iteration);
 *  - optional quasi-definite mode (Options::delta > 0): -delta on exactly-zero
 *    diagonal entries (multiplier block of a cg2o KKT system), so every
 *    diagonal pivot is acceptable and the symmetric pivot tolerance can be 0;
 *    adaptive refinement against the unregularized matrix removes the bias.
 *    Refinement runs only when some entry was regularized.
 * Use it when pivoting is needed (indefinite matrices that are not
 * quasi-definite). For SPD or quasi-definite systems LinearSolverEigenQDLDLT
 * is several times faster. Marginal covariances are not supported (LU).
 */
template <typename MatrixType>
class LinearSolverUmfpackFast : public g2o::LinearSolverCCS<MatrixType> {
public:
  struct Options {
    int strategy = UMFPACK_STRATEGY_SYMMETRIC; // or _AUTO, _UNSYMMETRIC
    int ordering = UMFPACK_ORDERING_AMD;
    double pivotTolerance = 0.1;      // UMFPACK default 0.1
    double symPivotTolerance = 0.001; // UMFPACK default 0.001
    int scale = UMFPACK_SCALE_SUM;    // UMFPACK default
    int umfpackRefinement = 0;        // UMFPACK_IRSTEP (UMFPACK default 2)
    double delta = 0.0;               // > 0 enables quasi-definite mode
    int refinementSteps =
        5; // max own refinement steps (only if entries were regularized)
    double refinementTol =
        1e-12; // stop when max|r| <= tol * max|b| (or stalled)
  };

  LinearSolverUmfpackFast() : g2o::LinearSolverCCS<MatrixType>() {
    umfpack_di_defaults(_control);
  }
  ~LinearSolverUmfpackFast() override { release(); }
  LinearSolverUmfpackFast(const LinearSolverUmfpackFast &) = delete;
  LinearSolverUmfpackFast &operator=(const LinearSolverUmfpackFast &) = delete;

  Options &options() { return _opt; }
  //! diagonal entries that received -delta in the last factorization
  int lastRegularizedCount() const { return _numRegularized; }

  bool init() override {
    _init = true; // pattern is re-checked; symbolic kept if unchanged
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
    const double t0 = g2o::get_monotonic_time();
    if (!fillMatrix(A))
      return false;
    applyControl();
    if (_symbolic == nullptr) {
      int status = umfpack_di_symbolic(_n, _n, _Ap.data(), _Ai.data(),
                                       _Ax.data(), &_symbolic, _control, _info);
      if (status != UMFPACK_OK) {
        G2O_ERROR("umfpack_di_symbolic failed ({})", status);
        return false;
      }
      setStat(true, g2o::get_monotonic_time() - t0);
    }
    const double t1 = g2o::get_monotonic_time();
    if (_numeric)
      umfpack_di_free_numeric(&_numeric);
    int status = umfpack_di_numeric(_Ap.data(), _Ai.data(), _Ax.data(),
                                    _symbolic, &_numeric, _control, _info);
    if (status != UMFPACK_OK) {
      G2O_ERROR("umfpack_di_numeric failed ({})", status);
      return false;
    }
    status = umfpack_di_solve(UMFPACK_A, _Ap.data(), _Ai.data(), _Ax.data(), x,
                              b, _numeric, _control, _info);
    if (status != UMFPACK_OK) {
      G2O_ERROR("umfpack_di_solve failed ({})", status);
      return false;
    }
    if (_numRegularized > 0)
      refine(x, b); // nothing to undo otherwise
    setStat(false, g2o::get_monotonic_time() - t1);
    for (int i = 0; i < _n; ++i)
      if (!std::isfinite(x[i])) {
        G2O_ERROR("LinearSolverUmfpackFast: solve produced non-finite values");
        return false;
      }
    return true;
  }

  Options _opt;
  bool _init = true;
  int _n = 0;
  double _control[UMFPACK_CONTROL];
  double _info[UMFPACK_INFO];
  void *_symbolic = nullptr;
  void *_numeric = nullptr;
  // g2o's upper-triangular CCS (source) and the full CSC used by UMFPACK
  std::vector<int> _Cp, _Ci, _Ap, _Ai, _mapLo, _mapUp, _diagPos;
  std::vector<double> _Cx, _Ax, _r, _dx;
  std::vector<char> _regMask; // diagonal entries that received -delta
  int _numRegularized = 0;    // how many of them in the last fill

  void release() {
    if (_numeric)
      umfpack_di_free_numeric(&_numeric);
    if (_symbolic)
      umfpack_di_free_symbolic(&_symbolic);
  }

  void applyControl() {
    _control[UMFPACK_STRATEGY] = _opt.strategy;
    _control[UMFPACK_ORDERING] = _opt.ordering;
    _control[UMFPACK_PIVOT_TOLERANCE] = _opt.pivotTolerance;
    _control[UMFPACK_SYM_PIVOT_TOLERANCE] = _opt.symPivotTolerance;
    _control[UMFPACK_SCALE] = _opt.scale;
    _control[UMFPACK_IRSTEP] = _opt.umfpackRefinement;
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

  // Build (once per pattern) the full CSC pattern and the scatter maps.
  bool fillMatrix(const g2o::SparseBlockMatrix<MatrixType> &A) {
    if (_init) {
      this->initMatrixStructure(A);
      const int n = A.rows(), nz = A.nonZeros();
      std::vector<int> Cp(n + 1), Ci(nz);
      _Cx.resize(nz);
      // A.nonZeros() counts whole blocks; the upper triangle has fewer entries
      const int nzU =
          this->_ccsMatrix->fillCCS(Cp.data(), Ci.data(), _Cx.data(), true);
      Ci.resize(nzU);
      _Cx.resize(nzU);
      const bool same =
          (_symbolic != nullptr) && n == _n && Cp == _Cp && Ci == _Ci;
      _init = false;
      if (!same) {
        release();
        _n = n;
        _Cp.swap(Cp);
        _Ci.swap(Ci);
        buildFullPattern();
      }
    } else {
      this->_ccsMatrix->fillCCS(_Cx.data(), true);
    }
    // scatter values
    const int nzU = static_cast<int>(_Cx.size());
    for (int e = 0; e < nzU; ++e) {
      _Ax[_mapUp[e]] = _Cx[e];
      if (_mapLo[e] >= 0)
        _Ax[_mapLo[e]] = _Cx[e];
    }
    _numRegularized = 0;
    if (_opt.delta > 0.0)
      for (int c = 0; c < _n; ++c) {
        _regMask[c] = _diagPos[c] >= 0 && _Ax[_diagPos[c]] == 0.0;
        if (_regMask[c]) {
          _Ax[_diagPos[c]] = -_opt.delta;
          ++_numRegularized;
        }
      }
    return true;
  }

  void buildFullPattern() {
    const int n = _n, nzU = static_cast<int>(_Ci.size());
    std::vector<int> cnt(n, 0);
    for (int j = 0; j < n; ++j)
      for (int k = _Cp[j]; k < _Cp[j + 1]; ++k) {
        ++cnt[j];
        if (_Ci[k] != j)
          ++cnt[_Ci[k]];
      }
    _Ap.assign(n + 1, 0);
    for (int j = 0; j < n; ++j)
      _Ap[j + 1] = _Ap[j] + cnt[j];
    // collect (row, source entry, upper?) per column, then sort rows
    std::vector<std::vector<std::pair<int, int>>> col(n); // (row, code)
    for (int j = 0; j < n; ++j)
      for (int k = _Cp[j]; k < _Cp[j + 1]; ++k) {
        const int i = _Ci[k];
        col[j].push_back({i, 2 * k}); // upper entry (i, j) in column j
        if (i != j)
          col[i].push_back({j, 2 * k + 1}); // mirrored (j, i) in column i
      }
    _Ai.resize(_Ap[n]);
    _Ax.assign(_Ap[n], 0.0);
    _mapUp.assign(nzU, -1);
    _mapLo.assign(nzU, -1);
    _diagPos.assign(n, -1);
    for (int j = 0; j < n; ++j) {
      std::sort(col[j].begin(), col[j].end());
      int p = _Ap[j];
      for (auto &rc : col[j]) {
        _Ai[p] = rc.first;
        const int k = rc.second / 2;
        if (rc.second % 2 == 0)
          _mapUp[k] = p;
        else
          _mapLo[k] = p;
        if (rc.first == j)
          _diagPos[j] = p;
        ++p;
      }
    }
    _r.resize(n);
    _dx.resize(n);
    _regMask.assign(n, 0);
  }

  // adaptive refinement against the unregularised matrix (see QDLDLT wrapper)
  void refine(double *x, const double *b) {
    double normB = 0.0, prev = 1e300;
    for (int i = 0; i < _n; ++i)
      normB = std::max(normB, std::abs(b[i]));
    for (int it = 0; it < _opt.refinementSteps; ++it) {
      for (int i = 0; i < _n; ++i)
        _r[i] = b[i];
      for (int j = 0; j < _n; ++j)
        for (int k = _Ap[j]; k < _Ap[j + 1]; ++k) {
          const double v =
              (k == _diagPos[j] && _regMask[j]) ? 0.0 : _Ax[k]; // undo -delta
          _r[_Ai[k]] -= v * x[j];
        }
      double rn = 0.0;
      for (int i = 0; i < _n; ++i)
        rn = std::max(rn, std::abs(_r[i]));
      if (rn <= _opt.refinementTol * normB || rn > 0.5 * prev)
        break;
      prev = rn;
      umfpack_di_solve(UMFPACK_A, _Ap.data(), _Ai.data(), _Ax.data(),
                       _dx.data(), _r.data(), _numeric, _control, _info);
      for (int i = 0; i < _n; ++i)
        x[i] += _dx[i];
    }
  }

  bool solveBlocks_impl(
      const g2o::SparseBlockMatrix<MatrixType> &,
      std::function<void(g2o::MarginalCovarianceCholesky &)>) override {
    G2O_WARN("LinearSolverUmfpackFast does not support marginal covariances");
    return false;
  }
};

} // namespace cg2o

#endif // CG2O_LINEAR_SOLVER_UMFPACK_FAST_H
