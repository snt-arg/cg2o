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

#ifndef G2O_LINEAR_SOLVER_EIGEN_LSCG_H
#define G2O_LINEAR_SOLVER_EIGEN_LSCG_H

#include <Eigen/IterativeLinearSolvers> // Include LeastSquaresConjugateGradient
#include <Eigen/Sparse>
#include <cassert>

#include "g2o/core/batch_stats.h"
#include "g2o/core/linear_solver.h"
#include "g2o/core/marginal_covariance_cholesky.h"
#include "g2o/stuff/logger.h"
#include "g2o/stuff/timeutil.h"

namespace cg2o {

/**
 * \brief Linear solver which uses the LeastSquaresConjugateGradient solver from
 * Eigen
 *
 * Implements an iterative solver using Eigen's LeastSquaresConjugateGradient
 * method.
 */
template <typename MatrixType>
class LinearSolverEigenLSCG : public g2o::LinearSolverCCS<MatrixType> {
public:
  typedef Eigen::SparseMatrix<double, Eigen::ColMajor> SparseMatrix;
  typedef Eigen::LeastSquaresConjugateGradient<SparseMatrix> LSCGDecomposition;

public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  LinearSolverEigenLSCG() : g2o::LinearSolverCCS<MatrixType>(), _init(true) {}

  virtual bool init() override {
    _init = true;
    return true;
  }

  virtual bool solve(const g2o::SparseBlockMatrix<MatrixType> &A, double *x,
                     double *b) override {
    using namespace g2o;

    double t;
    if (!computeLSCG(A, t))
      return false;

    // Solving the linear system using LeastSquaresConjugateGradient
    VectorX::MapType xx(x, _sparseMatrix.cols());
    VectorX::ConstMapType bb(b, _sparseMatrix.cols());
    xx = _lscg.solve(bb);

    if (_lscg.info() != Eigen::Success) {
      G2O_ERROR("LeastSquaresConjugateGradient solve failed.");
      return false;
    }

    // Collect statistics
    G2OBatchStatistics *globalStats = G2OBatchStatistics::globalStats();
    if (globalStats) {
      globalStats->timeNumericDecomposition = get_monotonic_time() - t;
    }

    return true;
  }

protected:
  bool _init;
  SparseMatrix _sparseMatrix;
  SparseMatrix _symmetricMatrix;
  LSCGDecomposition _lscg;

  // Compute LeastSquaresConjugateGradient decomposition
  bool computeLSCG(const g2o::SparseBlockMatrix<MatrixType> &A, double &t) {
    if (_init)
      _sparseMatrix.resize(A.rows(), A.cols());
    fillSparseMatrix(A, !_init);
    _init = false;

    t = g2o::get_monotonic_time();
    _symmetricMatrix =
        _sparseMatrix.selfadjointView<Eigen::Upper>(); // Enforce symmetry
    _lscg.compute(_symmetricMatrix);

    if (_lscg.info() != Eigen::Success) {
      G2O_ERROR("LeastSquaresConjugateGradient decomposition failed.");
      return false;
    }

    return true;
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
    double t;
    if (!computeLSCG(A, t))
      return false;
    return true;
  }
};

} // namespace cg2o

#endif
