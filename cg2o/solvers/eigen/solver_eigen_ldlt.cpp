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
#include "linear_solver_eigen_ldlt.h"

using namespace std;

namespace cg2o {

namespace {
template <int p, int l, bool blockorder>
 std::unique_ptr<g2o::BlockSolverBase> AllocateSolverLDLT() {
  G2O_DEBUG(
      "Using EigenLDLT poseDim {} landMarkDim {} blockordering {}", p,
      l, blockorder);
  auto linearSolver = std::make_unique<
      LinearSolverEigenLDLT<typename g2o::BlockSolverPL<p, l>::PoseMatrixType>>();
  linearSolver->setBlockOrdering(blockorder);
  return std::make_unique<g2o::BlockSolverPL<p, l>>(std::move(linearSolver));
}
}  // namespace

/**
 * Helper function for allocating LDLT-based solvers
 */
static g2o::OptimizationAlgorithm* createSolverLDLT(const std::string& ldltSolverName) {
  static const std::map<std::string,
                        std::function< std::unique_ptr<g2o::BlockSolverBase>()>>
      solver_factories{
          {"var", &AllocateSolverLDLT<-1, -1, true>},
          {"fix3_2", &AllocateSolverLDLT<3, 2, true>},
          {"fix6_3", &AllocateSolverLDLT<6, 3, true>},
          {"fix7_3", &AllocateSolverLDLT<7, 3, true>},
          {"fix3_2_scalar", &AllocateSolverLDLT<3, 2, false>},
          {"fix6_3_scalar", &AllocateSolverLDLT<6, 3, false>},
          {"fix7_3_scalar", &AllocateSolverLDLT<7, 3, false>},
      };
  using namespace g2o;
  string solverName = ldltSolverName.substr(3);
  auto solverf = solver_factories.find(solverName);
  if (solverf == solver_factories.end()) return nullptr;

  string methodName = ldltSolverName.substr(0, 2);

  if (methodName == "gn") {
    return new OptimizationAlgorithmGaussNewton(solverf->second());
  } else if (methodName == "lm") {
    return new OptimizationAlgorithmLevenberg(solverf->second());
  } else if (methodName == "dl") {
    return new OptimizationAlgorithmDogleg(solverf->second());
  }

  return nullptr;
}

class EigenLDLT_SolverCreator : public g2o::AbstractOptimizationAlgorithmCreator {
 public:
  explicit EigenLDLT_SolverCreator(const g2o::OptimizationAlgorithmProperty& p)
      : AbstractOptimizationAlgorithmCreator(p) {}
  virtual g2o::OptimizationAlgorithm* construct() {
    return createSolverLDLT(property().name);
  }
};

// Register LDLT-based solvers
G2O_REGISTER_OPTIMIZATION_LIBRARY(eigen_ldlt);

G2O_REGISTER_OPTIMIZATION_ALGORITHM(gn_var_ldlt, new EigenLDLT_SolverCreator(
    g2o::OptimizationAlgorithmProperty("gn_var_ldlt", "Gauss-Newton: LDLT solver with Eigen's SimplicialLDLT (variable blocksize)", "Eigen", false, Eigen::Dynamic, Eigen::Dynamic)));

G2O_REGISTER_OPTIMIZATION_ALGORITHM(lm_var_ldlt, new EigenLDLT_SolverCreator(
    g2o::OptimizationAlgorithmProperty("lm_var_ldlt", "Levenberg: LDLT solver with Eigen's SimplicialLDLT (variable blocksize)", "Eigen", false, Eigen::Dynamic, Eigen::Dynamic)));

}  // namespace cg2o
