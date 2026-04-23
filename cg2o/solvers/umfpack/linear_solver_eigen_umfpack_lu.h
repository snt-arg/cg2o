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

#ifndef G2O_LINEAR_SOLVER_EIGEN_UMFPACK_LU_H
#define G2O_LINEAR_SOLVER_EIGEN_UMFPACK_LU_H

#include <Eigen/Sparse>
#include <Eigen/UmfPackSupport>
#include <Eigen/src/Core/util/Constants.h>
#include <cassert>

#include "g2o/core/batch_stats.h"
#include "g2o/core/linear_solver.h"
#include "g2o/core/marginal_covariance_cholesky.h"
#include "g2o/stuff/logger.h"
#include "g2o/stuff/timeutil.h"

#include <filesystem>

namespace cg2o {

/**
 * \brief Linear solver using UmfPackLU (SuiteSparse) from Eigen
 *
 * Implements a robust linear solver using Eigen's UmfPackLU decomposition.
 */
template <typename MatrixType>
class LinearSolverEigenUmfPackLU : public g2o::LinearSolverCCS<MatrixType> {
public:
  typedef Eigen::SparseMatrix<double, Eigen::ColMajor> SparseMatrix;
  typedef Eigen::Triplet<double> Triplet;

  using UmfPackLUDecomposition = Eigen::UmfPackLU<SparseMatrix>;

public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  LinearSolverEigenUmfPackLU()
      : g2o::LinearSolverCCS<MatrixType>(), _init(true) {}

  virtual bool init() override {
    _init = true;
    return true;
  }

  void debug_linear_solver(const auto &A, const double *b = nullptr,
                           const double *x = nullptr) {
    using namespace g2o;
    namespace fs = std::filesystem;

    std::string folder_name = "umfpack_linear_debug";

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

  virtual bool solve(const g2o::SparseBlockMatrix<MatrixType> &A, double *x,
                     double *b) override {
    using namespace g2o;
    double t;
    // auto t1 = get_monotonic_time();
    if (!computeUmfPackLU(A, t))
      return false;

    // Solve linear system using UmfPackLU
    VectorX::MapType xx(x, _sparseMatrix.cols());
    VectorX::ConstMapType bb(b, _sparseMatrix.cols());
    xx = _umfpackLU.solve(bb);
    // auto t2 = get_monotonic_time();
    // std::cout << "UmfPackLU solve time: " << t2 - t1 << std::endl;

    if (_umfpackLU.info() != Eigen::Success) {
      G2O_ERROR("UmfPackLU solve failed.");
      debug_linear_solver(A, b, x);
      return false;
    }

#if CG2O_DEBUG_LINEAR_SOLVER
    debug_linear_solver(A, b, x);
#endif

    // Collect statistics
    G2OBatchStatistics *globalStats = G2OBatchStatistics::globalStats();
    if (globalStats) {
      globalStats->timeNumericDecomposition = get_monotonic_time() - t;
      globalStats->choleskyNNZ = 0;
    }

    return true;
  }

protected:
  bool _init;
  SparseMatrix _sparseMatrix;
  SparseMatrix _symmetricMatrix;
  UmfPackLUDecomposition _umfpackLU;

  // Compute UmfPackLU decomposition
  bool computeUmfPackLU(const g2o::SparseBlockMatrix<MatrixType> &A,
                        double &t) {
    using namespace g2o;
    // Resize matrix if needed
    if (_init)
      _sparseMatrix.resize(A.rows(), A.cols());
    // auto t1 = get_monotonic_time();
    fillSparseMatrix(A, !_init);
    // auto t2 = get_monotonic_time();
    // std::cout << "Fill sparse matrix time: " << t2 - t1 << std::endl;
    if (_init)
      computeSymbolicDecomposition();
    _init = false;

    t = get_monotonic_time();

    // Compute LU decomposition with selfadjointView
    _symmetricMatrix = _sparseMatrix.selfadjointView<Eigen::Upper>();
    _umfpackLU.factorize(_symmetricMatrix);

    if (_umfpackLU.info() != Eigen::Success) {
      G2O_ERROR("UmfPackLU decomposition failed.");
      debug_linear_solver(A, nullptr, nullptr);

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
    double t = get_monotonic_time();
    _symmetricMatrix = _sparseMatrix.selfadjointView<Eigen::Upper>();
    _umfpackLU.analyzePattern(_symmetricMatrix);
    _init = false; // Mark as initialized

    G2OBatchStatistics *globalStats = G2OBatchStatistics::globalStats();
    if (globalStats)
      globalStats->timeSymbolicDecomposition = get_monotonic_time() - t;
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
    if (!computeUmfPackLU(A, t))
      return false;

    g2o::MarginalCovarianceCholesky mcc;
    return true;
  }
};

} // namespace cg2o

#endif // G2O_LINEAR_SOLVER_EIGEN_UMFPACK_LU_H
