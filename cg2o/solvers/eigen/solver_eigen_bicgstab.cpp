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

#include "g2o/core/block_solver.h"
#include "g2o/core/optimization_algorithm_dogleg.h"
#include "g2o/core/optimization_algorithm_factory.h"
#include "g2o/core/optimization_algorithm_gauss_newton.h"
#include "g2o/core/optimization_algorithm_levenberg.h"
#include "g2o/core/sparse_optimizer.h"
#include "g2o/stuff/logger.h"
#include "linear_solver_eigen_bicgstab.h" // Updated header

using namespace std;

namespace cg2o {

namespace {
template <int p, int l, bool blockorder>
std::unique_ptr<g2o::BlockSolverBase> AllocateSolverBiCGSTAB() {
  G2O_DEBUG("Using EigenBiCGSTAB poseDim {} landMarkDim {} blockordering {}", p,
            l, blockorder);
  auto linearSolver = std::make_unique<LinearSolverEigenBiCGSTAB<
      typename g2o::BlockSolverPL<p, l>::PoseMatrixType>>();
  linearSolver->setBlockOrdering(blockorder);
  return std::make_unique<g2o::BlockSolverPL<p, l>>(std::move(linearSolver));
}
} // namespace

/**
 * Helper function for allocating BiCGSTAB-based solvers
 */
static g2o::OptimizationAlgorithm *
createSolverBiCGSTAB(const std::string &fullSolverName) {

  static const std::map<std::string,
                        std::function<std::unique_ptr<g2o::BlockSolverBase>()>>
      solver_factories{
          {"var", &AllocateSolverBiCGSTAB<-1, -1, true>},
          {"fix3_2", &AllocateSolverBiCGSTAB<3, 2, true>},
          {"fix6_3", &AllocateSolverBiCGSTAB<6, 3, true>},
          {"fix7_3", &AllocateSolverBiCGSTAB<7, 3, true>},
          {"fix3_2_scalar", &AllocateSolverBiCGSTAB<3, 2, false>},
          {"fix6_3_scalar", &AllocateSolverBiCGSTAB<6, 3, false>},
          {"fix7_3_scalar", &AllocateSolverBiCGSTAB<7, 3, false>},
      };
  using namespace g2o;
  string solverName = fullSolverName.substr(3);
  auto solverf = solver_factories.find(solverName);
  if (solverf == solver_factories.end())
    return nullptr;

  string methodName = fullSolverName.substr(0, 2);

  if (methodName == "gn") {
    return new OptimizationAlgorithmGaussNewton(solverf->second());
  } else if (methodName == "lm") {
    return new OptimizationAlgorithmLevenberg(solverf->second());
  } else if (methodName == "dl") {
    return new OptimizationAlgorithmDogleg(solverf->second());
  }

  return nullptr;
}

class EigenBiCGSTAB_SolverCreator
    : public g2o::AbstractOptimizationAlgorithmCreator {
public:
  explicit EigenBiCGSTAB_SolverCreator(
      const g2o::OptimizationAlgorithmProperty &p)
      : AbstractOptimizationAlgorithmCreator(p) {}
  virtual g2o::OptimizationAlgorithm *construct() {
    return createSolverBiCGSTAB(property().name);
  }
};

// Register BiCGSTAB-based solvers
G2O_REGISTER_OPTIMIZATION_LIBRARY(eigen_bicgstab);

G2O_REGISTER_OPTIMIZATION_ALGORITHM(
    gn_var_bicgstab,
    new EigenBiCGSTAB_SolverCreator(g2o::OptimizationAlgorithmProperty(
        "gn_var_bicgstab",
        "Gauss-Newton: BiCGSTAB solver using Eigen's iterative methods "
        "(variable blocksize)",
        "Eigen", false, Eigen::Dynamic, Eigen::Dynamic)));
G2O_REGISTER_OPTIMIZATION_ALGORITHM(
    gn_fix3_2_bicgstab,
    new EigenBiCGSTAB_SolverCreator(g2o::OptimizationAlgorithmProperty(
        "gn_fix3_2_bicgstab",
        "Gauss-Newton: BiCGSTAB solver using Eigen's iterative methods (fixed "
        "blocksize)",
        "Eigen", true, 3, 2)));
G2O_REGISTER_OPTIMIZATION_ALGORITHM(
    gn_fix6_3_bicgstab,
    new EigenBiCGSTAB_SolverCreator(g2o::OptimizationAlgorithmProperty(
        "gn_fix6_3_bicgstab",
        "Gauss-Newton: BiCGSTAB solver using Eigen's iterative methods (fixed "
        "blocksize)",
        "Eigen", true, 6, 3)));
G2O_REGISTER_OPTIMIZATION_ALGORITHM(
    gn_fix7_3_bicgstab,
    new EigenBiCGSTAB_SolverCreator(g2o::OptimizationAlgorithmProperty(
        "gn_fix7_3_bicgstab",
        "Gauss-Newton: BiCGSTAB solver using Eigen's iterative methods (fixed "
        "blocksize)",
        "Eigen", true, 7, 3)));

G2O_REGISTER_OPTIMIZATION_ALGORITHM(
    lm_var_bicgstab,
    new EigenBiCGSTAB_SolverCreator(g2o::OptimizationAlgorithmProperty(
        "lm_var_bicgstab",
        "Levenberg: BiCGSTAB solver using Eigen's iterative methods (variable "
        "blocksize)",
        "Eigen", false, Eigen::Dynamic, Eigen::Dynamic)));
G2O_REGISTER_OPTIMIZATION_ALGORITHM(
    lm_fix3_2_bicgstab,
    new EigenBiCGSTAB_SolverCreator(g2o::OptimizationAlgorithmProperty(
        "lm_fix3_2_bicgstab",
        "Levenberg: BiCGSTAB solver using Eigen's iterative methods (fixed "
        "blocksize)",
        "Eigen", true, 3, 2)));
G2O_REGISTER_OPTIMIZATION_ALGORITHM(
    lm_fix6_3_bicgstab,
    new EigenBiCGSTAB_SolverCreator(g2o::OptimizationAlgorithmProperty(
        "lm_fix6_3_bicgstab",
        "Levenberg: BiCGSTAB solver using Eigen's iterative methods (fixed "
        "blocksize)",
        "Eigen", true, 6, 3)));
G2O_REGISTER_OPTIMIZATION_ALGORITHM(
    lm_fix7_3_bicgstab,
    new EigenBiCGSTAB_SolverCreator(g2o::OptimizationAlgorithmProperty(
        "lm_fix7_3_bicgstab",
        "Levenberg: BiCGSTAB solver using Eigen's iterative methods (fixed "
        "blocksize)",
        "Eigen", true, 7, 3)));

G2O_REGISTER_OPTIMIZATION_ALGORITHM(
    dl_var_bicgstab,
    new EigenBiCGSTAB_SolverCreator(g2o::OptimizationAlgorithmProperty(
        "dl_var_bicgstab",
        "Dogleg: BiCGSTAB solver using Eigen's iterative methods (variable "
        "blocksize)",
        "Eigen", false, Eigen::Dynamic, Eigen::Dynamic)));

} // namespace cg2o
