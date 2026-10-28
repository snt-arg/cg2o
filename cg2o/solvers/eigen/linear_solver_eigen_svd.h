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

#ifndef G2O_LINEAR_SOLVER_EIGEN_SVD_H
#define G2O_LINEAR_SOLVER_EIGEN_SVD_H

#include <Eigen/Core>
#include <Eigen/SVD>
#include <cassert>
#include <filesystem>
#include <macros.h>

#include "g2o/core/linear_solver.h"

#include <algorithm>

#include "linear_solver_checks.h"

namespace cg2o {

/**
 * \brief linear solver using dense SVD decomposition
 */
template <typename MatrixType>
class LinearSolverEigenSVD : public g2o::LinearSolver<MatrixType> {
public:
  LinearSolverEigenSVD() : g2o::LinearSolver<MatrixType>() {}

  virtual ~LinearSolverEigenSVD() {}

  virtual bool init() { return true; }

  void debug_linear_solver(const auto &A, const double *b = nullptr,
                           const double *x = nullptr) {
    using namespace g2o;
    namespace fs = std::filesystem;
    std::string folder_name = "svd_linear_debug";

    // 1. Create folder if it doesn't exist
    if (!fs::exists(folder_name)) {
      fs::create_directory(folder_name);
    }

    // 2. Find the next available index 'x_val'
    int x_val = 0;
    while (fs::exists(folder_name + "/" + std::to_string(x_val) + "_Hessian" +
                      ".txt")) {
      x_val++;
    }

    // 3. Define paths for the new files
    std::string base_path = folder_name + "/";
    std::string hessian_file =
        base_path + std::to_string(x_val) + "_Hessian" + ".txt";
    std::string rhs_file = base_path + std::to_string(x_val) + "_rhs" + ".txt";
    std::string sol_file =
        base_path + std::to_string(x_val) + "_solution" + ".txt";

    // 4. Save data (Assumes your existing write functions work)
    A.writeOctave(hessian_file.c_str());
    if (b) {
      writeVector(rhs_file.c_str(), b, A.cols());
    }
    if (x) {
      writeVector(sol_file.c_str(), x, A.cols());
    }

    std::cout << "Debug files saved to: " << folder_name << " (Index: " << x_val
              << ")" << std::endl;
  }

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

    // Solving via SVD decomposition
    _decomposition.compute(_H, Eigen::ComputeThinU | Eigen::ComputeThinV);
    VectorX::MapType xvec(x, n);
    VectorX::ConstMapType bvec(b, n);

    xvec = _decomposition.solve(bvec);

#if CG2O_DEBUG_LINEAR_SOLVER
    debug_linear_solver(A, b, x);
#endif

    return _decomposition.info() == Eigen::Success &&
           internal::acceptableDenseSolution(_H, xvec, bvec);
  }

protected:
  g2o::MatrixX _H;
  Eigen::BDCSVD<g2o::MatrixX> _decomposition;
};

} // namespace cg2o

#endif
