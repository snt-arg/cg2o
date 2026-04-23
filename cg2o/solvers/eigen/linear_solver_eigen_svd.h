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

namespace cg2o {

/**
 * \brief linear solver using dense SVD decomposition
 */
template <typename MatrixType>
class LinearSolverEigenSVD : public g2o::LinearSolver<MatrixType> {
public:
  LinearSolverEigenSVD() : g2o::LinearSolver<MatrixType>(), _reset(true) {}

  virtual ~LinearSolverEigenSVD() {}

  virtual bool init() {
    _reset = true;
    return true;
  }

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
    using namespace g2o;

    int n = A.cols();
    int m = A.cols();

    MatrixX &H = _H;
    if (H.cols() != n) {
      H.resize(n, m);
      _reset = true;
    }
    if (_reset) {
      _reset = false;
      H.setZero();
    }

    // Copy the sparse block matrix into a dense matrix
    int c_idx = 0;
    for (size_t i = 0; i < A.blockCols().size(); ++i) {
      int c_size = A.colsOfBlock(i);
      assert(c_idx == A.colBaseOfBlock(i) && "Mismatch in block indices");

      const typename SparseBlockMatrix<MatrixType>::IntBlockMap &col =
          A.blockCols()[i];
      if (col.size() > 0) {
        typename SparseBlockMatrix<MatrixType>::IntBlockMap::const_iterator it;
        for (it = col.begin(); it != col.end(); ++it) {
          int r_idx = A.rowBaseOfBlock(it->first);
          // Only the upper triangular block is processed
          if (it->first <= (int)i) {
            int r_size = A.rowsOfBlock(it->first);
            H.block(r_idx, c_idx, r_size, c_size) = *(it->second);
            if (r_idx != c_idx) // Write the lower triangular block
              H.block(c_idx, r_idx, c_size, r_size) = it->second->transpose();
          }
        }
      }

      c_idx += c_size;
    }

    // Solving via SVD decomposition
    Eigen::BDCSVD<MatrixX> svd(H, Eigen::ComputeThinU | Eigen::ComputeThinV);
    VectorX::MapType xvec(x, m);
    VectorX::ConstMapType bvec(b, n);

    xvec = svd.solve(bvec);

#if CG2O_DEBUG_LINEAR_SOLVER
    debug_linear_solver(A, b, x);
#endif

    return svd.info() == Eigen::Success;
  }

protected:
  bool _reset;
  g2o::MatrixX _H;
};

} // namespace cg2o

#endif
