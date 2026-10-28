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

#ifndef G2O_LINEAR_SOLVER_EIGEN_DENSE_PARTIAL_LU_H
#define G2O_LINEAR_SOLVER_EIGEN_DENSE_PARTIAL_LU_H

#include <Eigen/Core>
#include <Eigen/LU>
#include <cassert>
#include <utility>
#include <vector>

#include "g2o/core/batch_stats.h"
#include "g2o/core/linear_solver.h"

#include <algorithm>

#include "linear_solver_checks.h"

namespace cg2o {

/**
 * \brief linear solver using dense Partial LU decomposition
 */
template <typename MatrixType>
class LinearSolverEigenDensePartialLU : public g2o::LinearSolver<MatrixType> {
public:
  LinearSolverEigenDensePartialLU() : g2o::LinearSolver<MatrixType>() {}

  virtual ~LinearSolverEigenDensePartialLU() {}

  virtual bool init() { return true; }

  bool solve(const g2o::SparseBlockMatrix<MatrixType> &A, double *x,
             double *b) {
    if (A.rows() == 0)
      return true;
    if (A.rows() == A.cols() && solveImpl(A, x, b))
      return true;
    std::fill(x, x + A.cols(), 0.0);
    return false;
  }

protected:
  bool solveImpl(const g2o::SparseBlockMatrix<MatrixType> &A, double *x,
                 double *b) {
    using namespace g2o;
    const int n = A.cols();
    internal::fillSymmetricDenseMatrix(A, _H);

    // Solving via Partial LU decomposition
    _decomposition.compute(_H);
    VectorX::MapType xvec(x, n);
    VectorX::ConstMapType bvec(b, n);

    xvec = _decomposition.solve(bvec);
    return internal::acceptableDenseSolution(_H, xvec, bvec);
  }

protected:
  g2o::MatrixX _H;
  Eigen::PartialPivLU<g2o::MatrixX> _decomposition;
};

} // namespace cg2o

#endif
