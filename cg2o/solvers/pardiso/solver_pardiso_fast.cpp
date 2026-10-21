/*
Contributed to cg2o; distributed under the same BSD 3-Clause license as the
rest of the library.

Registers LinearSolverPardisoFast with g2o's OptimizationAlgorithmFactory, so it can be
created by name: <method>_<blocksize>_pardiso_fast, e.g. "gn_var_pardiso_fast"
(method gn | lm | dl; blocksize var | fix3_2 | fix6_3 | fix7_3, with
"_scalar" for scalar instead of block ordering). With a static cg2o build,
put G2O_USE_OPTIMIZATION_LIBRARY(pardiso_fast) into the program that uses the
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
#include "linear_solver_pardiso_fast.h"

namespace cg2o {

namespace {
template <int p, int l, bool blockorder>
std::unique_ptr<g2o::BlockSolverBase> AllocateSolverPardisoFast() {
  G2O_DEBUG("Using LinearSolverPardisoFast poseDim {} landMarkDim {} blockordering {}", p, l,
            blockorder);
  auto linearSolver = std::make_unique<
      LinearSolverPardisoFast<typename g2o::BlockSolverPL<p, l>::PoseMatrixType>>();
  linearSolver->setBlockOrdering(blockorder);
  return std::make_unique<g2o::BlockSolverPL<p, l>>(std::move(linearSolver));
}

g2o::OptimizationAlgorithm* createSolverPardisoFast(const std::string& fullSolverName) {
  static const std::map<std::string,
                        std::function<std::unique_ptr<g2o::BlockSolverBase>()>>
      solver_factories{
          {"var", &AllocateSolverPardisoFast<-1, -1, true>},
          {"fix3_2", &AllocateSolverPardisoFast<3, 2, true>},
          {"fix6_3", &AllocateSolverPardisoFast<6, 3, true>},
          {"fix7_3", &AllocateSolverPardisoFast<7, 3, true>},
          {"var_scalar", &AllocateSolverPardisoFast<-1, -1, false>},
          {"fix3_2_scalar", &AllocateSolverPardisoFast<3, 2, false>},
          {"fix6_3_scalar", &AllocateSolverPardisoFast<6, 3, false>},
          {"fix7_3_scalar", &AllocateSolverPardisoFast<7, 3, false>},
      };
  // <method>_<blocksize>_pardiso_fast
  const std::string suffix = "_pardiso_fast";
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

class PardisoFast_SolverCreator : public g2o::AbstractOptimizationAlgorithmCreator {
 public:
  explicit PardisoFast_SolverCreator(const g2o::OptimizationAlgorithmProperty& p)
      : g2o::AbstractOptimizationAlgorithmCreator(p) {}
  g2o::OptimizationAlgorithm* construct() override {
    return createSolverPardisoFast(property().name);
  }
};
}  // namespace

G2O_REGISTER_OPTIMIZATION_LIBRARY(pardiso_fast);

#define CG2O_REGISTER_PardisoFast(method, blocks, desc, fixed, p, l)            \
  G2O_REGISTER_OPTIMIZATION_ALGORITHM(                                     \
      method##_##blocks##_pardiso_fast,                                        \
      new PardisoFast_SolverCreator(g2o::OptimizationAlgorithmProperty(          \
          #method "_" #blocks "_pardiso_fast", desc, "PARDISO", fixed, p, l)))

CG2O_REGISTER_PardisoFast(gn, var, "Gauss-Newton: MKL PARDISO symmetric LDL^T with cached analysis (variable blocksize)", false, Eigen::Dynamic, Eigen::Dynamic);
CG2O_REGISTER_PardisoFast(gn, fix3_2, "Gauss-Newton: MKL PARDISO symmetric LDL^T with cached analysis (fixed blocksize)", true, 3, 2);
CG2O_REGISTER_PardisoFast(gn, fix6_3, "Gauss-Newton: MKL PARDISO symmetric LDL^T with cached analysis (fixed blocksize)", true, 6, 3);
CG2O_REGISTER_PardisoFast(gn, fix7_3, "Gauss-Newton: MKL PARDISO symmetric LDL^T with cached analysis (fixed blocksize)", true, 7, 3);
CG2O_REGISTER_PardisoFast(lm, var, "Levenberg: MKL PARDISO symmetric LDL^T with cached analysis (variable blocksize)", false, Eigen::Dynamic, Eigen::Dynamic);
CG2O_REGISTER_PardisoFast(lm, fix3_2, "Levenberg: MKL PARDISO symmetric LDL^T with cached analysis (fixed blocksize)", true, 3, 2);
CG2O_REGISTER_PardisoFast(lm, fix6_3, "Levenberg: MKL PARDISO symmetric LDL^T with cached analysis (fixed blocksize)", true, 6, 3);
CG2O_REGISTER_PardisoFast(lm, fix7_3, "Levenberg: MKL PARDISO symmetric LDL^T with cached analysis (fixed blocksize)", true, 7, 3);
CG2O_REGISTER_PardisoFast(dl, var, "Dogleg: MKL PARDISO symmetric LDL^T with cached analysis (variable blocksize)", false, Eigen::Dynamic, Eigen::Dynamic);
CG2O_REGISTER_PardisoFast(dl, fix3_2, "Dogleg: MKL PARDISO symmetric LDL^T with cached analysis (fixed blocksize)", true, 3, 2);
CG2O_REGISTER_PardisoFast(dl, fix6_3, "Dogleg: MKL PARDISO symmetric LDL^T with cached analysis (fixed blocksize)", true, 6, 3);
CG2O_REGISTER_PardisoFast(dl, fix7_3, "Dogleg: MKL PARDISO symmetric LDL^T with cached analysis (fixed blocksize)", true, 7, 3);

#undef CG2O_REGISTER_PardisoFast

}  // namespace cg2o
