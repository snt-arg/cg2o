/*
Copyright (c) 2023, University of Luxembourg
All rights reserved.

Redistributions and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its
   contributors may be used to endorse or promote products derived from
   this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS 'AS IS'
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
*/

#ifndef G2O_LINEAR_SOLVER_EIGEN_LDLT_H
#define G2O_LINEAR_SOLVER_EIGEN_LDLT_H

#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>
#include <Eigen/src/Core/util/Constants.h>
#include <algorithm>
#include <cassert>
#include <vector>

#include "g2o/core/batch_stats.h"
#include "g2o/core/linear_solver.h"
#include "g2o/core/marginal_covariance_cholesky.h"
#include "g2o/stuff/logger.h"
#include "g2o/stuff/timeutil.h"

namespace cg2o {

/**
 * \brief Linear solver which uses the sparse LDLT solver from Eigen
 *
 * Implements a robust linear solver using Eigen's SimplicialLDLT decomposition.
 */
template <typename MatrixType>
class LinearSolverEigenLDLT : public g2o::LinearSolverCCS<MatrixType> {
public:
  typedef Eigen::SparseMatrix<double, Eigen::ColMajor> SparseMatrix;
  typedef Eigen::Triplet<double> Triplet;

  // Upper: g2o fills only the upper triangle, so Eigen can use it directly.
  // (With the default Lower, passing selfadjointView<Upper>() made Eigen build
  // a temporary full copy of the matrix on every factorization.)
  using ldltDecomposition =
      Eigen::SimplicialLDLT<SparseMatrix, Eigen::Upper, Eigen::AMDOrdering<int>>;

public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  LinearSolverEigenLDLT() : g2o::LinearSolverCCS<MatrixType>(), _init(true) {}

  virtual bool init() override {
    _init = true;
    return true;
  }

  virtual bool solve(const g2o::SparseBlockMatrix<MatrixType> &A, double *x,
                     double *b) override {
    using namespace g2o;
    double t;
    // g2o's GN/LM apply the update vector even when solve() fails, so a
    // failed solve leaves x = 0 (no step) instead of stale values.
    if (!computeLDLT(A, t)) {
      std::fill(x, x + A.rows(), 0.0);
      return false;
    }

    // Solving the linear system using LDLT decomposition
    VectorX::MapType xx(x, _sparseMatrix.cols());
    VectorX::ConstMapType bb(b, _sparseMatrix.cols());
    xx = _ldlt.solve(bb);

    if (_ldlt.info() != Eigen::Success) {
      G2O_ERROR("LDLT solve failed.");
      xx.setZero();
      return false;
    }

    // Collect statistics
    G2OBatchStatistics *globalStats = G2OBatchStatistics::globalStats();
    if (globalStats) {
      globalStats->timeNumericDecomposition = get_monotonic_time() - t;
      globalStats->choleskyNNZ = _ldlt.matrixL().nestedExpression().nonZeros();
    }

    return true;
  }

protected:
  bool _init;
  SparseMatrix _sparseMatrix;
  ldltDecomposition _ldlt;
  // sparsity pattern of the last symbolic analysis (see computeLDLT)
  bool _analyzed = false;
  std::vector<int> _patternOuter, _patternInner;

  // Compute LDLT decomposition
  bool computeLDLT(const g2o::SparseBlockMatrix<MatrixType> &A, double &t) {
    // Perform some operations only once upon init
    if (_init)
      _sparseMatrix.resize(A.rows(), A.cols());
    fillSparseMatrix(A, !_init);
    // g2o calls init() at every optimize(). Redo the ordering and symbolic
    // analysis only if the sparsity pattern actually changed (in MPC it never
    // does), instead of at every optimize().
    if (_init && !samePattern())
      computeSymbolicDecomposition();
    _init = false;

    t = g2o::get_monotonic_time();
    _ldlt.factorize(_sparseMatrix); // upper triangle, see ldltDecomposition
    if (_ldlt.info() != Eigen::Success) {
      if (this->writeDebug()) {
        G2O_ERROR("LDLT decomposition failed, writing debug.txt");
        A.writeOctave("debug.txt");
      } else {
        G2O_DEBUG("LDLT decomposition failed");
      }
      return false;
    }

    return true;
  }

  /**
   * Compute the symbolic decomposition of the matrix only once.
   * Since A has the same pattern in all the iterations, we only
   * compute the fill-in reducing ordering once and re-use for all
   * the following iterations.
   */
  void computeSymbolicDecomposition() {
    using namespace g2o;
    double t = g2o::get_monotonic_time();
    _ldlt.analyzePattern(_sparseMatrix);
    storePattern();
    _init = false; // Mark as initialized

    G2OBatchStatistics *globalStats = G2OBatchStatistics::globalStats();
    if (globalStats)
      globalStats->timeSymbolicDecomposition = get_monotonic_time() - t;
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

  void fillSparseMatrix(const g2o::SparseBlockMatrix<MatrixType> &A,
                        bool onlyValues) {
    if (onlyValues) {
      this->_ccsMatrix->fillCCS(_sparseMatrix.valuePtr(), true);
      return;
    }
    this->initMatrixStructure(A);
    _sparseMatrix.resizeNonZeros(A.nonZeros());
    int nz = this->_ccsMatrix->fillCCS(_sparseMatrix.outerIndexPtr(),
                                       _sparseMatrix.innerIndexPtr(),
                                       _sparseMatrix.valuePtr(), true);
    (void)nz;
    assert(nz <= static_cast<int>(_sparseMatrix.data().size()));
  }

  virtual bool solveBlocks_impl(
      const g2o::SparseBlockMatrix<MatrixType> &A,
      std::function<void(g2o::MarginalCovarianceCholesky &)> /*compute*/)
      override {
    // Marginal covariances are not implemented here. Returning true (as
    // before) made SparseOptimizer::computeMarginals report success with an
    // empty result. Use LinearSolverEigenQDLDLT, which supports them for SPD
    // systems.
    (void)A;
    G2O_WARN("LinearSolverEigenLDLT does not support marginal covariances; "
             "use LinearSolverEigenQDLDLT");
    return false;
  }
};

} // namespace cg2o

#endif
