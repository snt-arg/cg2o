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

#ifndef G2O_LINEAR_SOLVER_EIGEN_DENSE_QR_H
#define G2O_LINEAR_SOLVER_EIGEN_DENSE_QR_H

#include <Eigen/Core>
#include <Eigen/QR>
#include <cassert>

#include "g2o/core/linear_solver.h"

#include <algorithm>

#include "linear_solver_checks.h"

namespace cg2o {

/**
 * \brief Linear solver using Eigen's dense QR decomposition
 */
template <typename MatrixType>
class LinearSolverEigenDenseQR : public g2o::LinearSolver<MatrixType> {
public:
  LinearSolverEigenDenseQR() : g2o::LinearSolver<MatrixType>() {}

  virtual ~LinearSolverEigenDenseQR() {}

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

    // Solving via QR decomposition
    _decomposition.compute(_H);
    VectorX::MapType xvec(x, n);
    VectorX::ConstMapType bvec(b, n);

    xvec = _decomposition.solve(bvec);
    // HouseholderQR has no info(); validate the answer against the system.
    return internal::acceptableDenseSolution(_H, xvec, bvec);
  }

protected:
  g2o::MatrixX _H;
  Eigen::HouseholderQR<g2o::MatrixX> _decomposition;
};

} // namespace cg2o

#endif // G2O_LINEAR_SOLVER_EIGEN_DENSE_QR_H
