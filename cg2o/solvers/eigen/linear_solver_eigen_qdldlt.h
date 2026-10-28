/*
Contributed to cg2o; distributed under the same BSD 3-Clause license as the
rest of the library (see the header of linear_solver_eigen_ldlt.h).
*/
#ifndef CG2O_LINEAR_SOLVER_EIGEN_QDLDLT_H
#define CG2O_LINEAR_SOLVER_EIGEN_QDLDLT_H

#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>
#include <Eigen/SparseLU>
#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <stdexcept>
#include <vector>

#include "g2o/core/batch_stats.h"
#include "g2o/core/linear_solver.h"
#include "g2o/core/marginal_covariance_cholesky.h"
#include "g2o/stuff/logger.h"
#include "g2o/stuff/timeutil.h"

namespace cg2o {

/**
 * Sparse LDL^T with optional diagonal regularization, iterative refinement,
 * and a residual-verified SparseLU fallback. Header-only; Eigen and g2o only.
 *
 * The incoming matrix K is stored as its upper triangle. Exactly-zero stored
 * diagonal entries receive -delta, forming K_reg = K + R. This assumes that
 * such entries belong to dual variables; it does not make an arbitrary
 * indefinite matrix quasi-definite. In a KKT system, quasi-definiteness also
 * requires a positive-definite primal block and negative-definite dual block.
 * Any regularization already applied by the caller remains part of K.
 *
 * Accuracy is measured against K, removing only this solver's own R:
 *     ||b - K*x||_inf <= fallbackAbsoluteTolerance()
 *                       + fallbackTolerance() * ||b||_inf.
 * Default: relative tolerance 1e-8, absolute tolerance 0. Setting the relative
 * tolerance to 0 explicitly disables verification AND the LU fallback. In
 * that mode, finite LDL^T iterates are accepted after optional refinement.
 *
 * Refinement corrections are enabled after internal regularization, while
 * verification is enabled (fallbackTolerance() > 0, the default), or by
 * setAlwaysRefine(true). The same current residual serves refinement and
 * verification. Correction iterates are accepted only if their residual is
 * finite and strictly decreases. The final step is never reported successful
 * with an excessive residual when verification is enabled.
 *
 * LDL^T symbolic analysis and the pre-permuted value map are reused while the
 * upper-triangle pattern stays fixed. The full symmetric LU pattern and its
 * upper-to-full value map are built lazily on the first fallback and reused.
 * LU numeric factorization is repeated whenever fallback is needed: changing
 * matrix values must not reuse old numeric factors. No values are pruned.
 *
 * Call init() whenever the sparsity pattern changes, including changes with
 * the same size and number of entries. Dimension changes are also detected
 * automatically. Instances are stateful and not thread-safe.
 *
 * A failed solve zeros x, because some g2o optimizers can apply x even after
 * a false return. x and b must be distinct nonoverlapping arrays. SPD marginal
 * covariances continue to use LDL^T and do not use the LU fallback.
 */
template <typename MatrixType>
class LinearSolverEigenQDLDLT : public g2o::LinearSolverCCS<MatrixType> {
public:
  using SparseMatrix = Eigen::SparseMatrix<double, Eigen::ColMajor, int>;
  using DecompositionBase = Eigen::SimplicialLDLT<SparseMatrix, Eigen::Upper,
                                                  Eigen::AMDOrdering<int>>;
  using VectorMap = Eigen::Map<Eigen::VectorXd>;
  using ConstVectorMap = Eigen::Map<const Eigen::VectorXd>;
  /**
   * SimplicialLDLT whose numeric factorization starts from a cached,
   * already permuted copy of the matrix.
   *
   * Eigen's factorize() rebuilds P A P^T (a sparse symmetric permutation:
   * pattern work, a column count and an allocation) at EVERY call although the
   * sparsity pattern is fixed between calls. In the closed-loop NMPC profile of
   * the px4 slalom this copy was 13 % of the linear solve (7.7 % of all
   * instructions of a solve).
   *
   * Cache the sparsity pattern of P A P^T and copy only numeric entries
   * before factorize_preordered(). This preserves Eigen's symbolic ordering.
   * This optional fast path uses protected Eigen internals: test it with the
   * Eigen version used by the application. setPrePermutedFactorization(false)
   * selects the public Eigen factorize() path instead.
   */
  class Decomposition : public DecompositionBase {
  public:
    //! to be called after analyzePattern(A); A must keep its pattern afterwards
    void buildPermutationMap(const SparseMatrix &A) {
      _haveMap = false;
      if (this->m_P.size() == 0)
        return; // no ordering: stock path
      SparseMatrix idx = A;
      for (Eigen::Index k = 0; k < idx.nonZeros(); ++k)
        idx.valuePtr()[k] = static_cast<double>(k + 1);
      _permuted.resize(A.rows(), A.cols());
      _permuted.template selfadjointView<Eigen::Upper>() =
          idx.template selfadjointView<Eigen::Upper>().twistedBy(this->m_P);
      _map.resize(_permuted.nonZeros());
      for (Eigen::Index k = 0; k < _permuted.nonZeros(); ++k)
        _map[k] = static_cast<int>(_permuted.valuePtr()[k]) - 1;
      _haveMap = true;
    }
    //! numeric factorization of A, which has the pattern given to
    //! buildPermutationMap()
    void factorizeFast(const SparseMatrix &A) {
      if (!_haveMap ||
          _permuted.nonZeros() != static_cast<Eigen::Index>(_map.size())) {
        this->factorize(A);
        return;
      }
      const double *v = A.valuePtr();
      double *pv = _permuted.valuePtr();
      const int *m = _map.data();
      const Eigen::Index nz = _permuted.nonZeros();
      for (Eigen::Index k = 0; k < nz; ++k)
        pv[k] = v[m[k]];
#if EIGEN_VERSION_AT_LEAST(3, 4, 90)
      // Eigen master: factorize_preordered<DoLDLT, NonHermitian>
      this->template factorize_preordered<true, false>(_permuted);
#else
      this->template factorize_preordered<true>(_permuted);
#endif
    }

  private:
    SparseMatrix _permuted;
    std::vector<int> _map;
    bool _haveMap = false;
  };

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  LinearSolverEigenQDLDLT() : g2o::LinearSolverCCS<MatrixType>() {}
  LinearSolverEigenQDLDLT(double delta, bool refine, int steps)
      : LinearSolverEigenQDLDLT() {
    setDualRegularization(delta);
    setAlwaysRefine(refine);
    setRefinementSteps(steps);
  }
  LinearSolverEigenQDLDLT(double delta, bool refine)
      : LinearSolverEigenQDLDLT(delta, refine, 5) {}

  bool init() override {
    _init = true;
    return true;
  }

  void setDualRegularization(double delta) {
    requireNonnegativeFinite(delta, "dual regularization");
    _delta = delta;
  }
  double dualRegularization() const { return _delta; }
  //! maximum iterative-refinement steps (default 5). One fixed step is not
  //! always enough: on ill-conditioned KKT matrices (large barrier weights)
  //! each step removes only a factor delta / (lambda_min + delta) of the bias,
  //! lambda_min being the smallest eigenvalue of J H^-1 J^T.
  void setRefinementSteps(int steps) {
    if (steps < 0)
      throw std::invalid_argument("negative refinement steps");
    _refinementSteps = steps;
  }
  int refinementSteps() const { return _refinementSteps; }
  void setRefinementTolerance(double tol) {
    requireNonnegativeFinite(tol, "refinement tolerance");
    _refinementTol = tol;
  }
  double refinementTolerance() const { return _refinementTol; }
  void setAlwaysRefine(bool always) { _alwaysRefine = always; }
  bool alwaysRefine() const { return _alwaysRefine; }
  void setPrePermutedFactorization(bool on) { _prePermuted = on; }
  bool prePermutedFactorization() const { return _prePermuted; }
  void setFallbackTolerance(double tol) {
    requireNonnegativeFinite(tol, "fallback tolerance");
    _fallbackTol = tol;
  }
  double fallbackTolerance() const { return _fallbackTol; }
  //! false (default): when neither LDL^T nor the LU reaches the tolerance, the
  //! answer with the smaller residual is returned as the step and
  //! inaccurateCount() is incremented; true: such a solve fails (x = 0, returns
  //! false)
  void setFailOnInaccurate(bool fail) { _failOnInaccurate = fail; }
  bool failOnInaccurate() const { return _failOnInaccurate; }
  void setFallbackAbsoluteTolerance(double tol) {
    requireNonnegativeFinite(tol, "fallback absolute tolerance");
    _fallbackAbsTol = tol;
  }
  double fallbackAbsoluteTolerance() const { return _fallbackAbsTol; }

  //! Attempted corrections, including a correction rejected for no improvement.
  int lastRefinementSteps() const { return _lastRefinementSteps; }
  int lastRegularizedCount() const { return _numRegularized; }
  //! All solve calls; solves whose answer came from the LU fallback (verified,
  //! or kept because it beat the LDL^T answer); LU attempts; verified LU
  //! successes; LU attempts that did not reach the tolerance.
  long solveCount() const { return _solveCount; }
  long fallbackCount() const { return _fallbackUsedCount; }
  long fallbackAttemptCount() const { return _fallbackAttemptCount; }
  long fallbackSuccessCount() const { return _fallbackSuccessCount; }
  long fallbackFailureCount() const {
    return _fallbackAttemptCount - _fallbackSuccessCount;
  }
  long failedSolveCount() const { return _failedSolveCount; }
  //! Calls rejected for excessive/nonfinite residual, including failed LU
  //! recovery.
  long inaccurateCount() const { return _inaccurateCount; }
  bool lastUsedFallback() const { return _lastUsedFallback; }
  //! Best candidate residual before failure zeroing; NaN if it was not
  //! measured.
  double lastResidualNorm() const { return _lastResidualNorm; }
  double lastRelativeResidual() const { return _lastRelativeResidual; }
  long residualEvaluationCount() const { return _residualEvaluations; }
  long luAnalysisCount() const { return _luAnalysisCount; }
  long luPatternBuildCount() const { return _luPatternBuildCount; }

  bool solve(const g2o::SparseBlockMatrix<MatrixType> &A, double *x,
             double *b) override {
    ++_solveCount;
    _lastRefinementSteps = 0;
    _lastUsedFallback = false;
    _lastResidualNorm = _lastRelativeResidual =
        std::numeric_limits<double>::quiet_NaN();
    const int n = A.rows();
    if (n != A.cols() || n < 0 || (n > 0 && (!x || !b || x == b)))
      return fail(x, n);
    if (n == 0) {
      _lastResidualNorm = _lastRelativeResidual = 0.0;
      return true;
    }
    VectorMap xx(x, n);
    ConstVectorMap bb(b, n);
    if (!bb.allFinite())
      return fail(x, n);
    const double normB = bb.cwiseAbs().maxCoeff();
    const bool ldltOK = factorize(A);
    const double t = g2o::get_monotonic_time();
    bool success = false;
    if (ldltOK) {
      xx = _ldlt.solve(bb);
      if (_fallbackTol > 0.0 || _numRegularized > 0 || _alwaysRefine)
        success = refineAndCheck(xx, bb, normB);
      else
        success = xx.allFinite();
    } else if (_matrixReady && _fallbackTol > 0.0) {
      // Assembly succeeded but numeric LDL^T failed. No LDL^T solve/refinement.
      recordResidual(std::numeric_limits<double>::infinity(), normB);
      success = solveWithLU(xx, bb, normB);
      if (!success) {
        ++_inaccurateCount;
        success = !_failOnInaccurate && std::isfinite(_lastResidualNorm) &&
                  xx.allFinite() &&
                  _lastResidualNorm < std::numeric_limits<double>::infinity();
      }
    }
    g2o::G2OBatchStatistics *stats = g2o::G2OBatchStatistics::globalStats();
    if (stats) {
      stats->timeNumericDecomposition += g2o::get_monotonic_time() - t;
      if (success)
        stats->choleskyNNZ =
            _lastUsedFallback ? _lu.nnzL() + _lu.nnzU()
                              : _ldlt.matrixL().nestedExpression().nonZeros();
    }
    if (!success)
      return fail(x, n);
    return true;
  }

protected:
  bool _init = true;
  bool _analyzed = false;
  bool _matrixReady = false;
  bool _alwaysRefine = false;
  bool _prePermuted = true;
  bool _failOnInaccurate = false;
  double _delta = 1e-9;
  int _refinementSteps = 5;
  double _refinementTol = 1e-12;
  int _lastRefinementSteps = 0;
  int _numRegularized = 0;
  double _fallbackTol = 1e-8;
  double _fallbackAbsTol = 0.0;
  long _fallbackUsedCount = 0;
  long _fallbackAttemptCount = 0;
  long _fallbackSuccessCount = 0;
  long _failedSolveCount = 0;
  long _inaccurateCount = 0;
  long _solveCount = 0;
  long _residualEvaluations = 0;
  long _luAnalysisCount = 0;
  long _luPatternBuildCount = 0;
  bool _lastUsedFallback = false;
  double _lastResidualNorm = std::numeric_limits<double>::quiet_NaN();
  double _lastRelativeResidual = std::numeric_limits<double>::quiet_NaN();
  bool _luAnalyzed = false;
  bool _fullPatternReady = false;
  Eigen::SparseLU<SparseMatrix, Eigen::COLAMDOrdering<int>> _lu;
  SparseMatrix _fullMatrix;
  std::vector<int> _fullToUpper, _fullDiag;
  SparseMatrix _sparseMatrix;
  Decomposition _ldlt;
  std::vector<int> _diagIndex;
  Eigen::VectorXd _regularization, _residual, _candidate, _candidateResidual;
  std::vector<int> _patternOuter, _patternInner;

  static void requireNonnegativeFinite(double value, const char *name) {
    if (!std::isfinite(value) || value < 0.0)
      throw std::invalid_argument(name);
  }

  bool fail(double *x, int n) {
    if (x && n > 0)
      std::fill(x, x + n, 0.0);
    ++_failedSolveCount;
    return false;
  }

  void recordResidual(double rn, double normB) {
    _lastResidualNorm = rn;
    _lastRelativeResidual =
        normB > 0.0
            ? rn / normB
            : (rn == 0.0 ? 0.0 : std::numeric_limits<double>::infinity());
  }

  bool acceptable(double rn, double normB) const {
    // Long double prevents overflow of tol*normB from silently accepting inf.
    return std::isfinite(rn) && static_cast<long double>(rn) <=
                                    static_cast<long double>(_fallbackAbsTol) +
                                        static_cast<long double>(_fallbackTol) *
                                            static_cast<long double>(normB);
  }

  template <typename Derived>
  double residualNorm(const Eigen::MatrixBase<Derived> &candidate,
                      const ConstVectorMap &bb, Eigen::VectorXd &residual) {
    ++_residualEvaluations;
    if (!candidate.allFinite())
      return std::numeric_limits<double>::infinity();
    residual.noalias() =
        _sparseMatrix.template selfadjointView<Eigen::Upper>() * candidate;
    residual.array() -= _regularization.array() * candidate.array();
    residual = bb - residual;
    if (!residual.allFinite())
      return std::numeric_limits<double>::infinity();
    return residual.cwiseAbs().maxCoeff();
  }

  bool refineAndCheck(VectorMap &xx, const ConstVectorMap &bb, double normB) {
    double rn = residualNorm(xx, bb, _residual);
    // also refine unregularized (e.g. SPD) systems while verification is on;
    // an accurate answer stops the loop before any correction
    const bool allowRefinement =
        _numRegularized > 0 || _alwaysRefine || _fallbackTol > 0.0;
    if (allowRefinement && std::isfinite(rn)) {
      for (int k = 0; k < _refinementSteps; ++k) {
        if (static_cast<long double>(rn) <=
            static_cast<long double>(_refinementTol) * normB)
          break;
        _candidate = _ldlt.solve(_residual);
        _candidate += xx;
        ++_lastRefinementSteps;
        const double next = residualNorm(_candidate, bb, _candidateResidual);
        if (!std::isfinite(next) || !(next < rn))
          break;
        xx = _candidate;
        _residual.swap(_candidateResidual);
        rn = next;
      }
    }
    recordResidual(rn, normB);
    if (_fallbackTol == 0.0)
      return xx.allFinite();
    if (acceptable(rn, normB))
      return true;
    if (solveWithLU(xx, bb, normB, rn))
      return true;
    ++_inaccurateCount;
    return !_failOnInaccurate && xx.allFinite();
  }

  void buildFullPattern() {
    // Index tags preserve explicit zeros and avoid coeffRef()/insertion later.
    SparseMatrix tags = _sparseMatrix;
    for (Eigen::Index i = 0; i < tags.nonZeros(); ++i)
      tags.valuePtr()[i] = static_cast<double>(i + 1);
    _fullMatrix = tags.template selfadjointView<Eigen::Upper>();
    _fullMatrix.makeCompressed();
    _fullToUpper.resize(_fullMatrix.nonZeros());
    _fullDiag.assign(_fullMatrix.cols(), -1);
    for (int c = 0; c < _fullMatrix.cols(); ++c) {
      for (int i = _fullMatrix.outerIndexPtr()[c];
           i < _fullMatrix.outerIndexPtr()[c + 1]; ++i) {
        _fullToUpper[i] = static_cast<int>(_fullMatrix.valuePtr()[i]) - 1;
        if (_fullMatrix.innerIndexPtr()[i] == c)
          _fullDiag[c] = i;
      }
    }
    _fullPatternReady = true;
    _luAnalyzed = false;
    ++_luPatternBuildCount;
  }

  bool solveWithLU(
      VectorMap &xx, const ConstVectorMap &bb, double normB,
      double currentResidual = std::numeric_limits<double>::infinity()) {
    ++_fallbackAttemptCount;
    if (!_fullPatternReady)
      buildFullPattern();
    const double *source = _sparseMatrix.valuePtr();
    double *dest = _fullMatrix.valuePtr();
    for (Eigen::Index i = 0; i < _fullMatrix.nonZeros(); ++i)
      dest[i] = source[_fullToUpper[i]];
    for (int c = 0; c < static_cast<int>(_fullDiag.size()); ++c)
      if (_fullDiag[c] >= 0)
        dest[_fullDiag[c]] -= _regularization[c];
    if (!_luAnalyzed) {
      _lu.analyzePattern(_fullMatrix);
      _luAnalyzed = true;
      ++_luAnalysisCount;
    }
    _lu.factorize(_fullMatrix);
    if (_lu.info() != Eigen::Success)
      return false;
    _candidate = _lu.solve(bb);
    if (_lu.info() != Eigen::Success || !_candidate.allFinite())
      return false;
    const double rn = residualNorm(_candidate, bb, _candidateResidual);
    if (!acceptable(rn, normB)) {
      if (std::isfinite(rn) &&
          !(rn >= currentResidual)) { // not accurate enough, but better than
                                      // the LDL^T answer
        xx = _candidate;
        _residual.swap(_candidateResidual);
        recordResidual(rn, normB);
        ++_fallbackUsedCount;
      }
      return false;
    }
    xx = _candidate;
    _residual.swap(_candidateResidual);
    recordResidual(rn, normB);
    ++_fallbackUsedCount;
    ++_fallbackSuccessCount;
    _lastUsedFallback = true;
    return true;
  }

  bool factorize(const g2o::SparseBlockMatrix<MatrixType> &A) {
    _matrixReady = false;
    if (A.rows() != A.cols() || A.rows() <= 0)
      return false;
    g2o::G2OBatchStatistics *stats = g2o::G2OBatchStatistics::globalStats();
    if (A.rows() != _sparseMatrix.rows() || A.cols() != _sparseMatrix.cols())
      _init = true;
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
        _ldlt.buildPermutationMap(_sparseMatrix);
        storePattern();
        _luAnalyzed = false;
        _fullPatternReady = false;
      }
      _init = false;
      if (stats)
        stats->timeSymbolicDecomposition = g2o::get_monotonic_time() - t;
    } else {
      this->_ccsMatrix->fillCCS(_sparseMatrix.valuePtr(), true);
      regularize();
    }
    // Reject invalid coefficients before entering either numeric factorization.
    if (!Eigen::Map<const Eigen::VectorXd>(_sparseMatrix.valuePtr(),
                                           _sparseMatrix.nonZeros())
             .allFinite())
      return false;
    _matrixReady = true;
    const double t = g2o::get_monotonic_time();
    if (_prePermuted)
      _ldlt.factorizeFast(_sparseMatrix);
    else
      _ldlt.factorize(_sparseMatrix);
    if (stats)
      stats->timeNumericDecomposition = g2o::get_monotonic_time() - t;
    if (_ldlt.info() != Eigen::Success) {
      if (this->writeDebug())
        A.writeOctave("debug.txt");
      G2O_DEBUG("QDLDLT numeric factorization failed");
      return false;
    }
    return true;
  }
  bool samePattern() const {
    if (!_analyzed)
      return false;
    const int n = static_cast<int>(_sparseMatrix.cols());
    if (static_cast<int>(_patternOuter.size()) != n + 1)
      return false;
    const int *outer = _sparseMatrix.outerIndexPtr();
    if (!std::equal(outer, outer + n + 1, _patternOuter.begin()))
      return false;
    const int nnz = outer[n];
    if (static_cast<int>(_patternInner.size()) != nnz)
      return false;
    const int *inner = _sparseMatrix.innerIndexPtr();
    return std::equal(inner, inner + nnz, _patternInner.begin());
  }

  void storePattern() {
    const int n = static_cast<int>(_sparseMatrix.cols());
    const int *outer = _sparseMatrix.outerIndexPtr();
    _patternOuter.assign(outer, outer + n + 1);
    const int *inner = _sparseMatrix.innerIndexPtr();
    _patternInner.assign(inner, inner + outer[n]);
    _analyzed = true;
  }

  void locateDiagonal() {
    const int n = static_cast<int>(_sparseMatrix.cols());
    _diagIndex.assign(n, -1);
    const int *outer = _sparseMatrix.outerIndexPtr();
    const int *inner = _sparseMatrix.innerIndexPtr();
    for (int c = 0; c < n; ++c)
      for (int k = outer[c]; k < outer[c + 1]; ++k)
        if (inner[k] == c) {
          _diagIndex[c] = k;
          break;
        }
    _regularization = Eigen::VectorXd::Zero(n);
    _residual.resize(n);
    _candidate.resize(n);
    _candidateResidual.resize(n);
  }

  void regularize() {
    double *values = _sparseMatrix.valuePtr();
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
      const g2o::SparseBlockMatrix<MatrixType> &A,
      std::function<void(g2o::MarginalCovarianceCholesky &)> compute) override {
    if (!factorize(A))
      return false;
    const Eigen::VectorXd &d = _ldlt.vectorD();
    if (_numRegularized > 0 || !d.allFinite() || (d.array() <= 0.0).any()) {
      G2O_WARN("LinearSolverEigenQDLDLT: marginal covariances need a positive "
               "definite system");
      return false;
    }
    const SparseMatrix &L = _ldlt.matrixL().nestedExpression(); // strict lower
    const int n = static_cast<int>(L.cols());
    std::vector<int> Lp(n + 1), Li(L.nonZeros() + n);
    std::vector<double> Lx(L.nonZeros() + n);
    int p = 0;
    for (int j = 0; j < n; ++j) {
      const double s = std::sqrt(d[j]);
      Lp[j] = p;
      Li[p] = j;
      Lx[p++] = s;
      for (SparseMatrix::InnerIterator it(L, j); it; ++it) { // rows sorted
        Li[p] = static_cast<int>(it.row());
        Lx[p++] = it.value() * s;
      }
    }
    Lp[n] = p;
    // Eigen: P A P^T is factored and P.indices()[i] is the new position of i,
    // i.e. the inverse permutation that MarginalCovarianceCholesky expects
    std::vector<int> pinv(n);
    const auto &P = _ldlt.permutationP();
    for (int i = 0; i < n; ++i)
      pinv[i] = P.size() == n ? P.indices()[i] : i;
    g2o::MarginalCovarianceCholesky mcc;
    mcc.setCholeskyFactor(n, Lp.data(), Li.data(), Lx.data(), pinv.data());
    compute(mcc);
    g2o::G2OBatchStatistics *globalStats =
        g2o::G2OBatchStatistics::globalStats();
    if (globalStats)
      globalStats->choleskyNNZ = p;
    return true;
  }
};

} // namespace cg2o

#endif // CG2O_LINEAR_SOLVER_EIGEN_QDLDLT_H
