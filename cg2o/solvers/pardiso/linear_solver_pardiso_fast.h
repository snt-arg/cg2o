/*
Contributed to cg2o; distributed under the same BSD 3-Clause license as the
rest of the library (see the header of linear_solver_eigen_pardiso_lu.h).
*/

#ifndef CG2O_LINEAR_SOLVER_PARDISO_FAST_H
#define CG2O_LINEAR_SOLVER_PARDISO_FAST_H

#include <mkl_pardiso.h>
#include <mkl_types.h>

#include <algorithm>
#include <cmath>
#include <cstring>
#include <functional>
#include <vector>

#include "g2o/core/batch_stats.h"
#include "g2o/core/linear_solver.h"
#include "g2o/core/marginal_covariance_cholesky.h"
#include "g2o/stuff/logger.h"
#include "g2o/stuff/timeutil.h"

namespace cg2o {

/**
 * \brief MKL PARDISO for any g2o system, using the symmetric storage g2o
 * assembles (upper triangle).
 *
 * Differences to LinearSolverEigenPardisoLU:
 *  - uses PARDISO's symmetric modes instead of unsymmetric LU (mtype 11):
 *    mtype -2 (default) symmetric indefinite LDL^T with Bunch-Kaufman pivoting
 *    inside supernodes, for cg2o KKT systems; mtype 2 Cholesky for SPD
 *    systems (pose graphs, bundle adjustment, AL). Only the upper triangle is
 *    stored and factorized;
 *  - calls the MKL C API directly; the upper CSR pattern PARDISO needs is
 *    built once from g2o's upper CCS and values are scattered into it, instead
 *    of Eigen copying/transposing the matrix on every factorization;
 *  - the reordering/symbolic phase (11) is kept across optimize() calls while
 *    the sparsity pattern is unchanged; each solve runs phases 22 and 33;
 *  - ordering, pivoting, scaling/matching, perturbation and refinement are
 *    configurable (iparm), see Options;
 *  - optional quasi-definite mode (Options::delta > 0): -delta on exactly-zero
 *    diagonal entries and plain 1x1 pivots, with adaptive refinement against
 *    the unregularized matrix (only when some entry was regularized).
 * PARDISO is multithreaded through MKL (MKL_NUM_THREADS); it pays off for
 * large or denser systems. Marginal covariances are not supported.
 */
template <typename MatrixType>
class LinearSolverPardisoFast : public g2o::LinearSolverCCS<MatrixType> {
 public:
  struct Options {
    MKL_INT mtype = -2;        // -2 symmetric indefinite, 2 SPD (upper-triangle storage only)
    MKL_INT ordering = 0;      // iparm[1]: 0 minimum degree, 2 METIS, 3 parallel METIS
    MKL_INT pivoting = 1;      // iparm[20]: 0 1x1, 1 Bunch-Kaufman 1x1/2x2 (symmetric indefinite)
    MKL_INT scaling = 0;       // iparm[10]
    MKL_INT matching = 0;      // iparm[12]
    MKL_INT perturbation = 8;  // iparm[9]: tiny pivots perturbed by 10^-value
    MKL_INT refinement = 0;    // iparm[7]: PARDISO's own refinement steps (0 = none)
    double delta = 0.0;        // > 0 enables quasi-definite mode
    int refinementSteps = 5;   // max own refinement steps (only if entries were regularized)
    double refinementTol = 1e-12;  // stop when max|r| <= tol * max|b| (or stalled)
  };

  LinearSolverPardisoFast() : g2o::LinearSolverCCS<MatrixType>() {
    std::memset(_pt, 0, sizeof(_pt));
  }
  ~LinearSolverPardisoFast() override { release(); }

  Options& options() { return _opt; }

  bool init() override {
    _init = true;  // pattern is re-checked; symbolic phase kept if unchanged
    return true;
  }

  bool solve(const g2o::SparseBlockMatrix<MatrixType>& A, double* x,
             double* b) override {
    const double t0 = g2o::get_monotonic_time();
    fillMatrix(A);
    MKL_INT error = 0, nrhs = 1;
    if (!_analyzed) {
      setupIparm();
      MKL_INT phase = 11;
      pardiso(_pt, &_maxfct, &_mnum, &_opt.mtype, &phase, &_n, _a.data(),
              _ia.data(), _ja.data(), nullptr, &nrhs, _iparm, &_msglvl,
              nullptr, nullptr, &error);
      if (error != 0) {
        G2O_ERROR("PARDISO analysis failed ({})", error);
        return false;
      }
      _analyzed = true;
      setStat(true, g2o::get_monotonic_time() - t0);
    }
    const double t1 = g2o::get_monotonic_time();
    MKL_INT phase = 22;
    pardiso(_pt, &_maxfct, &_mnum, &_opt.mtype, &phase, &_n, _a.data(),
            _ia.data(), _ja.data(), nullptr, &nrhs, _iparm, &_msglvl, nullptr,
            nullptr, &error);
    if (error != 0) {
      G2O_ERROR("PARDISO factorization failed ({})", error);
      return false;
    }
    phase = 33;
    pardiso(_pt, &_maxfct, &_mnum, &_opt.mtype, &phase, &_n, _a.data(),
            _ia.data(), _ja.data(), nullptr, &nrhs, _iparm, &_msglvl, b, x,
            &error);
    if (error != 0) {
      G2O_ERROR("PARDISO solve failed ({})", error);
      return false;
    }
    if (_numRegularized > 0) refine(x, b);  // nothing to undo otherwise
    setStat(false, g2o::get_monotonic_time() - t1);
    return true;
  }

  //! diagonal entries that received -delta in the last factorization
  int lastRegularizedCount() const { return _numRegularized; }
  //! number of pivots PARDISO had to perturb in the last factorization
  int perturbedPivots() const { return static_cast<int>(_iparm[13]); }

 protected:
  Options _opt;
  bool _init = true, _analyzed = false;
  void* _pt[64];
  MKL_INT _iparm[64] = {0};
  MKL_INT _maxfct = 1, _mnum = 1, _msglvl = 0, _n = 0;
  std::vector<MKL_INT> _ia, _ja;           // upper CSR (row pointers, columns)
  std::vector<int> _Cp, _Ci, _map, _diagPos;
  std::vector<double> _Cx, _a, _r, _dx;
  std::vector<char> _regMask;
  int _numRegularized = 0;  // entries that received -delta in the last fill

  void release() {
    if (_analyzed) {
      MKL_INT phase = -1, error = 0, nrhs = 1;
      pardiso(_pt, &_maxfct, &_mnum, &_opt.mtype, &phase, &_n, nullptr,
              _ia.data(), _ja.data(), nullptr, &nrhs, _iparm, &_msglvl,
              nullptr, nullptr, &error);
      _analyzed = false;
    }
  }

  void setupIparm() {
    pardisoinit(_pt, &_opt.mtype, _iparm);  // MKL defaults for this mtype
    _iparm[0] = 1;                          // use the values below
    _iparm[1] = _opt.ordering;
    _iparm[7] = _opt.refinement;
    _iparm[9] = _opt.perturbation;
    _iparm[10] = _opt.scaling;
    _iparm[12] = _opt.matching;
    _iparm[20] = _opt.delta > 0.0 ? 0 : _opt.pivoting;
    _iparm[26] = 0;                         // no matrix checker
    _iparm[34] = 1;                         // zero-based indexing
  }

  void setStat(bool symbolic, double dt) {
    g2o::G2OBatchStatistics* s = g2o::G2OBatchStatistics::globalStats();
    if (!s) return;
    if (symbolic)
      s->timeSymbolicDecomposition = dt;
    else
      s->timeNumericDecomposition = dt;
  }

  void fillMatrix(const g2o::SparseBlockMatrix<MatrixType>& A) {
    if (_init) {
      this->initMatrixStructure(A);
      const int n = A.rows(), nz = A.nonZeros();
      std::vector<int> Cp(n + 1), Ci(nz);
      _Cx.resize(nz);
      // A.nonZeros() counts whole blocks; the upper triangle has fewer entries
      const int nzU = this->_ccsMatrix->fillCCS(Cp.data(), Ci.data(), _Cx.data(), true);
      Ci.resize(nzU);
      _Cx.resize(nzU);
      const bool same = _analyzed && n == _n && Cp == _Cp && Ci == _Ci;
      _init = false;
      if (!same) {
        release();
        _n = n;
        _Cp.swap(Cp);
        _Ci.swap(Ci);
        buildUpperCsr();
      }
    } else {
      this->_ccsMatrix->fillCCS(_Cx.data(), true);
    }
    const int nzU = static_cast<int>(_Cx.size());
    for (int e = 0; e < nzU; ++e) _a[_map[e]] = _Cx[e];
    _numRegularized = 0;
    if (_opt.delta > 0.0)
      for (int r = 0; r < _n; ++r) {
        _regMask[r] = _diagPos[r] >= 0 && _a[_diagPos[r]] == 0.0;
        if (_regMask[r]) {
          _a[_diagPos[r]] = -_opt.delta;
          ++_numRegularized;
        }
      }
  }

  // upper CCS entry (i, j), i <= j, becomes upper CSR row i, column j
  void buildUpperCsr() {
    const int n = _n, nzU = static_cast<int>(_Ci.size());
    std::vector<std::vector<std::pair<int, int>>> row(n);  // (col, source entry)
    for (int j = 0; j < n; ++j)
      for (int k = _Cp[j]; k < _Cp[j + 1]; ++k) row[_Ci[k]].push_back({j, k});
    _ia.assign(n + 1, 0);
    _ja.resize(nzU);
    _a.assign(nzU, 0.0);
    _map.assign(nzU, -1);
    _diagPos.assign(n, -1);
    int p = 0;
    for (int i = 0; i < n; ++i) {
      std::sort(row[i].begin(), row[i].end());
      _ia[i] = p;
      for (auto& ce : row[i]) {
        _ja[p] = ce.first;
        _map[ce.second] = p;
        if (ce.first == i) _diagPos[i] = p;
        ++p;
      }
    }
    _ia[n] = p;
    _r.resize(n);
    _dx.resize(n);
    _regMask.assign(n, 0);
  }

  // adaptive refinement against the unregularised matrix (see QDLDLT wrapper)
  void refine(double* x, const double* b) {
    MKL_INT phase = 33, error = 0, nrhs = 1;
    double normB = 0.0, prev = 1e300;
    for (int i = 0; i < _n; ++i) normB = std::max(normB, std::abs(b[i]));
    for (int it = 0; it < _opt.refinementSteps; ++it) {
      for (int i = 0; i < _n; ++i) _r[i] = b[i];
      for (int i = 0; i < _n; ++i)  // symmetric product from the upper CSR
        for (MKL_INT k = _ia[i]; k < _ia[i + 1]; ++k) {
          const int j = _ja[k];
          const double v = (j == i && _regMask[i]) ? 0.0 : _a[k];
          _r[i] -= v * x[j];
          if (j != i) _r[j] -= v * x[i];
        }
      double rn = 0.0;
      for (int i = 0; i < _n; ++i) rn = std::max(rn, std::abs(_r[i]));
      if (rn <= _opt.refinementTol * normB || rn > 0.5 * prev) break;
      prev = rn;
      pardiso(_pt, &_maxfct, &_mnum, &_opt.mtype, &phase, &_n, _a.data(),
              _ia.data(), _ja.data(), nullptr, &nrhs, _iparm, &_msglvl,
              _r.data(), _dx.data(), &error);
      for (int i = 0; i < _n; ++i) x[i] += _dx[i];
    }
  }

  bool solveBlocks_impl(
      const g2o::SparseBlockMatrix<MatrixType>&,
      std::function<void(g2o::MarginalCovarianceCholesky&)>) override {
    G2O_WARN("LinearSolverPardisoFast does not support marginal covariances");
    return false;
  }
};

}  // namespace cg2o

#endif  // CG2O_LINEAR_SOLVER_PARDISO_FAST_H
