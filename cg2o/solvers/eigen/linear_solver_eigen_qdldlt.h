/*
Contributed to cg2o; distributed under the same BSD 3-Clause license as the
rest of the library (see the header of linear_solver_eigen_ldlt.h).
*/

#ifndef CG2O_LINEAR_SOLVER_EIGEN_QDLDLT_H
#define CG2O_LINEAR_SOLVER_EIGEN_QDLDLT_H

#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>
#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <vector>

#include "g2o/core/batch_stats.h"
#include "g2o/core/linear_solver.h"
#include "g2o/core/marginal_covariance_cholesky.h"
#include "g2o/stuff/logger.h"
#include "g2o/stuff/timeutil.h"

namespace cg2o {

/**
 * \brief Sparse LDL^T solver (Eigen SimplicialLDLT, AMD ordering) for any g2o
 * problem, with optional quasi-definite regularization for KKT systems.
 *
 * Header-only; depends only on Eigen and g2o. Works with BlockSolverX and the
 * fixed-size block solvers, with or without the Schur complement.
 *
 * SPD systems (pose graphs, bundle adjustment, cg2o AL, unconstrained):
 *   no diagonal entry is zero, nothing is regularized and no refinement is
 *   done, so this is a plain simplicial LDL^T. Marginal covariances
 *   (SparseOptimizer::computeMarginals) are supported for SPD systems.
 *
 * KKT systems (cg2o ISPD/BIPM, where the equality multipliers are vertices):
 *
 *     [ H   J^T ]
 *     [ J   0   ]
 *
 *   is symmetric indefinite: Cholesky (g2o's Eigen/CSparse/CHOLMOD-supernodal
 *   solvers) cannot factor it and LU (UMFPACK, PARDISO-LU) ignores the
 *   symmetry. Two ways to use this solver:
 *   1. Let the solver regularize: every diagonal entry that is exactly zero
 *      gets -delta (default 1e-9). The matrix is then symmetric
 *      quasi-definite, which has an LDL^T factorization for any symmetric
 *      ordering (Vanderbei 1995). Adaptive iterative refinement against the
 *      unregularized matrix removes the bias of delta.
 *   2. Let the equality edges regularize (optimizer->setEqualityRegularization
 *      or edge->setDualRegularization): the matrix arrives as
 *      [H J^T; J -sigma I], nothing is zero, and the solver only factors it.
 *
 * The fill-reducing ordering and symbolic factorization are computed once and
 * kept across optimize() calls while the sparsity pattern is unchanged (g2o
 * calls init() at every optimize(); in MPC the pattern never changes).
 *
 * Assumption for mode 1: an exactly-zero diagonal entry belongs to a dual
 * (multiplier) variable. A primal variable without any cost term (zero
 * diagonal) would also receive -delta; give it a small prior instead.
 * No pivoting is done: for indefinite matrices that are not quasi-definite
 * use LinearSolverPardisoFast (Bunch-Kaufman) or LinearSolverUmfpackFast.
 */
template <typename MatrixType>
class LinearSolverEigenQDLDLT : public g2o::LinearSolverCCS<MatrixType> {
 public:
  using SparseMatrix = Eigen::SparseMatrix<double, Eigen::ColMajor, int>;
  using Decomposition =
      Eigen::SimplicialLDLT<SparseMatrix, Eigen::Upper, Eigen::AMDOrdering<int>>;

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  LinearSolverEigenQDLDLT() : g2o::LinearSolverCCS<MatrixType>() {}

  bool init() override {
    _init = true;  // pattern is re-checked; symbolic kept if unchanged
    return true;
  }

  //! value put on exactly-zero diagonal entries as -delta (default 1e-9);
  //! 0 disables the regularization
  void setDualRegularization(double delta) { _delta = delta; }
  double dualRegularization() const { return _delta; }
  //! maximum iterative-refinement steps (default 5). Refinement runs only when
  //! entries were regularized (or always, see setAlwaysRefine) and stops once
  //! max|b - K x| <= tol * max|b| or the residual stops shrinking. One fixed
  //! step is not always enough: on ill-conditioned KKT matrices (large barrier
  //! weights) each step removes only a factor delta / (lambda_min + delta) of
  //! the bias, lambda_min being the smallest eigenvalue of J H^-1 J^T.
  void setRefinementSteps(int steps) { _refinementSteps = steps; }
  int refinementSteps() const { return _refinementSteps; }
  //! stop refining once max|b - K x| <= tol * max|b| (default 1e-12)
  void setRefinementTolerance(double tol) { _refinementTol = tol; }
  double refinementTolerance() const { return _refinementTol; }
  //! also refine when nothing was regularized (removes rounding error of the
  //! factorization; costs at least one sparse product per solve). Default off.
  void setAlwaysRefine(bool always) { _alwaysRefine = always; }
  bool alwaysRefine() const { return _alwaysRefine; }
  //! refinement steps taken in the last solve
  int lastRefinementSteps() const { return _lastRefinementSteps; }
  //! diagonal entries that received -delta in the last factorization
  int lastRegularizedCount() const { return _numRegularized; }

  bool solve(const g2o::SparseBlockMatrix<MatrixType>& A, double* x,
             double* b) override {
    if (!factorize(A)) return false;
    const double t = g2o::get_monotonic_time();
    Eigen::Map<Eigen::VectorXd> xx(x, _sparseMatrix.cols());
    Eigen::Map<const Eigen::VectorXd> bb(b, _sparseMatrix.cols());
    xx = _ldlt.solve(bb);
    _lastRefinementSteps = 0;
    if (_numRegularized > 0 || _alwaysRefine) refine(xx, bb);
    if (!xx.allFinite()) {
      G2O_ERROR("QDLDLT solve produced non-finite values");
      return false;
    }
    g2o::G2OBatchStatistics* globalStats =
        g2o::G2OBatchStatistics::globalStats();
    if (globalStats) {
      globalStats->timeNumericDecomposition += g2o::get_monotonic_time() - t;
      globalStats->choleskyNNZ = _ldlt.matrixL().nestedExpression().nonZeros();
    }
    return true;
  }

 protected:
  bool _init = true;
  bool _analyzed = false;
  bool _alwaysRefine = false;
  double _delta = 1e-7;
  int _refinementSteps = 3;
  double _refinementTol = 1e-12;
  int _lastRefinementSteps = 0;
  int _numRegularized = 0;
  SparseMatrix _sparseMatrix;  // upper triangle, as filled by g2o
  Decomposition _ldlt;
  std::vector<int> _diagIndex;      // position of each diagonal entry
  Eigen::VectorXd _regularization;  // -delta where applied, else 0
  Eigen::VectorXd _residual;
  std::vector<int> _patternOuter, _patternInner;

  // adaptive iterative refinement against the unregularized matrix
  void refine(Eigen::Map<Eigen::VectorXd>& xx,
              const Eigen::Map<const Eigen::VectorXd>& bb) {
    const double normB = bb.cwiseAbs().maxCoeff();
    double prev = std::numeric_limits<double>::infinity();
    for (int k = 0; k < _refinementSteps; ++k) {
      // K x = (K_reg - R) x, R = diag(-delta on regularized entries)
      _residual = bb - (_sparseMatrix.selfadjointView<Eigen::Upper>() * xx -
                        _regularization.cwiseProduct(xx));
      const double rn = _residual.cwiseAbs().maxCoeff();
      // relative to |b| on purpose: a normwise backward error |K||x| would
      // hide the delta bias when K has huge barrier entries
      if (rn <= _refinementTol * normB) break;
      if (rn > 0.5 * prev) break;  // stalled: rounding level reached
      prev = rn;
      xx += _ldlt.solve(_residual);
      ++_lastRefinementSteps;
    }
  }

  bool factorize(const g2o::SparseBlockMatrix<MatrixType>& A) {
    g2o::G2OBatchStatistics* globalStats =
        g2o::G2OBatchStatistics::globalStats();
    if (_init) {
      const double t = g2o::get_monotonic_time();
      _sparseMatrix.resize(A.rows(), A.cols());
      this->initMatrixStructure(A);
      _sparseMatrix.resizeNonZeros(A.nonZeros());
      this->_ccsMatrix->fillCCS(_sparseMatrix.outerIndexPtr(),
                                _sparseMatrix.innerIndexPtr(),
                                _sparseMatrix.valuePtr(), true);
      locateDiagonal();
      regularize();
      if (!samePattern()) {
        _ldlt.analyzePattern(_sparseMatrix);
        storePattern();
      }
      _init = false;
      if (globalStats)
        globalStats->timeSymbolicDecomposition = g2o::get_monotonic_time() - t;
    } else {
      this->_ccsMatrix->fillCCS(_sparseMatrix.valuePtr(), true);
      regularize();
    }
    const double t = g2o::get_monotonic_time();
    _ldlt.factorize(_sparseMatrix);
    if (globalStats)
      globalStats->timeNumericDecomposition = g2o::get_monotonic_time() - t;
    if (_ldlt.info() != Eigen::Success) {
      if (this->writeDebug()) {
        G2O_ERROR("QDLDLT factorization failed, writing debug.txt");
        A.writeOctave("debug.txt");
      } else {
        G2O_DEBUG("QDLDLT factorization failed");
      }
      return false;
    }
    return true;
  }

  bool samePattern() const {
    if (!_analyzed) return false;
    const int n = static_cast<int>(_sparseMatrix.cols());
    if (static_cast<int>(_patternOuter.size()) != n + 1) return false;
    const int* outer = _sparseMatrix.outerIndexPtr();
    if (!std::equal(outer, outer + n + 1, _patternOuter.begin())) return false;
    const int nnz = outer[n];
    if (static_cast<int>(_patternInner.size()) != nnz) return false;
    const int* inner = _sparseMatrix.innerIndexPtr();
    return std::equal(inner, inner + nnz, _patternInner.begin());
  }

  void storePattern() {
    const int n = static_cast<int>(_sparseMatrix.cols());
    const int* outer = _sparseMatrix.outerIndexPtr();
    _patternOuter.assign(outer, outer + n + 1);
    const int* inner = _sparseMatrix.innerIndexPtr();
    _patternInner.assign(inner, inner + outer[n]);
    _analyzed = true;
  }

  void locateDiagonal() {
    const int n = static_cast<int>(_sparseMatrix.cols());
    _diagIndex.assign(n, -1);
    const int* outer = _sparseMatrix.outerIndexPtr();
    const int* inner = _sparseMatrix.innerIndexPtr();
    for (int c = 0; c < n; ++c)
      for (int k = outer[c]; k < outer[c + 1]; ++k)
        if (inner[k] == c) {
          _diagIndex[c] = k;
          break;
        }
    _regularization = Eigen::VectorXd::Zero(n);
    _residual.resize(n);
  }

  void regularize() {
    double* values = _sparseMatrix.valuePtr();
    _numRegularized = 0;
    for (int c = 0; c < static_cast<int>(_diagIndex.size()); ++c) {
      const int k = _diagIndex[c];
      if (_delta > 0.0 && k >= 0 && values[k] == 0.0) {
        values[k] = -_delta;
        _regularization[c] = -_delta;
        ++_numRegularized;
      } else {
        _regularization[c] = 0.0;
      }
    }
  }

  /**
   * Marginal covariances for SPD systems: L D L^T = P A P^T is turned into
   * the Cholesky factor L sqrt(D) that g2o's MarginalCovarianceCholesky
   * expects (lower CCS, diagonal first in each column, inverse permutation).
   * Not available for indefinite (KKT) systems.
   */
  bool solveBlocks_impl(
      const g2o::SparseBlockMatrix<MatrixType>& A,
      std::function<void(g2o::MarginalCovarianceCholesky&)> compute)
      override {
    if (!factorize(A)) return false;
    const Eigen::VectorXd& d = _ldlt.vectorD();
    if (_numRegularized > 0 || (d.array() <= 0.0).any()) {
      G2O_WARN(
          "LinearSolverEigenQDLDLT: marginal covariances need a positive "
          "definite system");
      return false;
    }
    const SparseMatrix& L = _ldlt.matrixL().nestedExpression();  // strict lower
    const int n = static_cast<int>(L.cols());
    std::vector<int> Lp(n + 1), Li(L.nonZeros() + n);
    std::vector<double> Lx(L.nonZeros() + n);
    int p = 0;
    for (int j = 0; j < n; ++j) {
      const double s = std::sqrt(d[j]);
      Lp[j] = p;
      Li[p] = j;
      Lx[p++] = s;
      for (SparseMatrix::InnerIterator it(L, j); it; ++it) {  // rows sorted
        Li[p] = static_cast<int>(it.row());
        Lx[p++] = it.value() * s;
      }
    }
    Lp[n] = p;
    // Eigen: P A P^T is factored and P.indices()[i] is the new position of i,
    // i.e. the inverse permutation that MarginalCovarianceCholesky expects
    std::vector<int> pinv(n);
    const auto& P = _ldlt.permutationP();
    for (int i = 0; i < n; ++i) pinv[i] = P.size() == n ? P.indices()[i] : i;
    g2o::MarginalCovarianceCholesky mcc;
    mcc.setCholeskyFactor(n, Lp.data(), Li.data(), Lx.data(), pinv.data());
    compute(mcc);
    g2o::G2OBatchStatistics* globalStats =
        g2o::G2OBatchStatistics::globalStats();
    if (globalStats) globalStats->choleskyNNZ = p;
    return true;
  }
};

}  // namespace cg2o

#endif  // CG2O_LINEAR_SOLVER_EIGEN_QDLDLT_H
