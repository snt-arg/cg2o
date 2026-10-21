/*
Contributed to cg2o; distributed under the same BSD 3-Clause license as the
rest of the library.

Registers LinearSolverEigenQDLDLT with g2o's OptimizationAlgorithmFactory, so it can be
created by name: <method>_<blocksize>_qdldlt, e.g. "gn_var_qdldlt"
(method gn | lm | dl; blocksize var | fix3_2 | fix6_3 | fix7_3, with
"_scalar" for scalar instead of block ordering). With a static cg2o build,
put G2O_USE_OPTIMIZATION_LIBRARY(eigen_qdldlt) into the program that uses the
factory so that this file is linked.
*/

#include <functional>
#include <map>
#include <memory>
#include <string>

#include "g2o/core/block_solver.h"
#include "g2o/core/optimization_algorithm_dogleg.h"
#include "g2o/core/optimization_algorithm_factory.h"
#include "g2o/core/optimization_algorithm_gauss_newton.h"
#include "g2o/core/optimization_algorithm_levenberg.h"
#include "g2o/stuff/logger.h"
#include "linear_solver_eigen_qdldlt.h"

namespace cg2o {

namespace {
template <int p, int l, bool blockorder>
std::unique_ptr<g2o::BlockSolverBase> AllocateSolverQDLDLT() {
  G2O_DEBUG("Using LinearSolverEigenQDLDLT poseDim {} landMarkDim {} blockordering {}", p, l,
            blockorder);
  auto linearSolver = std::make_unique<
      LinearSolverEigenQDLDLT<typename g2o::BlockSolverPL<p, l>::PoseMatrixType>>();
  linearSolver->setBlockOrdering(blockorder);
  return std::make_unique<g2o::BlockSolverPL<p, l>>(std::move(linearSolver));
}

g2o::OptimizationAlgorithm* createSolverQDLDLT(const std::string& fullSolverName) {
  static const std::map<std::string,
                        std::function<std::unique_ptr<g2o::BlockSolverBase>()>>
      solver_factories{
          {"var", &AllocateSolverQDLDLT<-1, -1, true>},
          {"fix3_2", &AllocateSolverQDLDLT<3, 2, true>},
          {"fix6_3", &AllocateSolverQDLDLT<6, 3, true>},
          {"fix7_3", &AllocateSolverQDLDLT<7, 3, true>},
          {"var_scalar", &AllocateSolverQDLDLT<-1, -1, false>},
          {"fix3_2_scalar", &AllocateSolverQDLDLT<3, 2, false>},
          {"fix6_3_scalar", &AllocateSolverQDLDLT<6, 3, false>},
          {"fix7_3_scalar", &AllocateSolverQDLDLT<7, 3, false>},
      };
  // <method>_<blocksize>_qdldlt
  const std::string suffix = "_qdldlt";
  if (fullSolverName.size() <= 3 + suffix.size() ||
      fullSolverName.compare(fullSolverName.size() - suffix.size(),
                             suffix.size(), suffix) != 0)
    return nullptr;
  const std::string method = fullSolverName.substr(0, 2);
  const std::string blocks =
      fullSolverName.substr(3, fullSolverName.size() - 3 - suffix.size());
  auto solverf = solver_factories.find(blocks);
  if (solverf == solver_factories.end()) return nullptr;
  if (method == "gn")
    return new g2o::OptimizationAlgorithmGaussNewton(solverf->second());
  if (method == "lm")
    return new g2o::OptimizationAlgorithmLevenberg(solverf->second());
  if (method == "dl")
    return new g2o::OptimizationAlgorithmDogleg(solverf->second());
  return nullptr;
}

class QDLDLT_SolverCreator : public g2o::AbstractOptimizationAlgorithmCreator {
 public:
  explicit QDLDLT_SolverCreator(const g2o::OptimizationAlgorithmProperty& p)
      : g2o::AbstractOptimizationAlgorithmCreator(p) {}
  g2o::OptimizationAlgorithm* construct() override {
    return createSolverQDLDLT(property().name);
  }
};
}  // namespace

G2O_REGISTER_OPTIMIZATION_LIBRARY(eigen_qdldlt);

#define CG2O_REGISTER_QDLDLT(method, blocks, desc, fixed, p, l)            \
  G2O_REGISTER_OPTIMIZATION_ALGORITHM(                                     \
      method##_##blocks##_qdldlt,                                        \
      new QDLDLT_SolverCreator(g2o::OptimizationAlgorithmProperty(          \
          #method "_" #blocks "_qdldlt", desc, "Eigen", fixed, p, l)))

CG2O_REGISTER_QDLDLT(gn, var, "Gauss-Newton: quasi-definite LDL^T with Eigen's SimplicialLDLT (variable blocksize)", false, Eigen::Dynamic, Eigen::Dynamic);
CG2O_REGISTER_QDLDLT(gn, fix3_2, "Gauss-Newton: quasi-definite LDL^T with Eigen's SimplicialLDLT (fixed blocksize)", true, 3, 2);
CG2O_REGISTER_QDLDLT(gn, fix6_3, "Gauss-Newton: quasi-definite LDL^T with Eigen's SimplicialLDLT (fixed blocksize)", true, 6, 3);
CG2O_REGISTER_QDLDLT(gn, fix7_3, "Gauss-Newton: quasi-definite LDL^T with Eigen's SimplicialLDLT (fixed blocksize)", true, 7, 3);
CG2O_REGISTER_QDLDLT(lm, var, "Levenberg: quasi-definite LDL^T with Eigen's SimplicialLDLT (variable blocksize)", false, Eigen::Dynamic, Eigen::Dynamic);
CG2O_REGISTER_QDLDLT(lm, fix3_2, "Levenberg: quasi-definite LDL^T with Eigen's SimplicialLDLT (fixed blocksize)", true, 3, 2);
CG2O_REGISTER_QDLDLT(lm, fix6_3, "Levenberg: quasi-definite LDL^T with Eigen's SimplicialLDLT (fixed blocksize)", true, 6, 3);
CG2O_REGISTER_QDLDLT(lm, fix7_3, "Levenberg: quasi-definite LDL^T with Eigen's SimplicialLDLT (fixed blocksize)", true, 7, 3);
CG2O_REGISTER_QDLDLT(dl, var, "Dogleg: quasi-definite LDL^T with Eigen's SimplicialLDLT (variable blocksize)", false, Eigen::Dynamic, Eigen::Dynamic);
CG2O_REGISTER_QDLDLT(dl, fix3_2, "Dogleg: quasi-definite LDL^T with Eigen's SimplicialLDLT (fixed blocksize)", true, 3, 2);
CG2O_REGISTER_QDLDLT(dl, fix6_3, "Dogleg: quasi-definite LDL^T with Eigen's SimplicialLDLT (fixed blocksize)", true, 6, 3);
CG2O_REGISTER_QDLDLT(dl, fix7_3, "Dogleg: quasi-definite LDL^T with Eigen's SimplicialLDLT (fixed blocksize)", true, 7, 3);

#undef CG2O_REGISTER_QDLDLT

}  // namespace cg2o
