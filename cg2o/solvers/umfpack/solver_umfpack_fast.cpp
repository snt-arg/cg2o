/*
Contributed to cg2o; distributed under the same BSD 3-Clause license as the
rest of the library.

Registers LinearSolverUmfpackFast with g2o's OptimizationAlgorithmFactory, so it
can be created by name: <method>_<blocksize>_umfpack_fast, e.g.
"gn_var_umfpack_fast" (method gn | lm | dl; blocksize var | fix3_2 | fix6_3 |
fix7_3, with
"_scalar" for scalar instead of block ordering). With a static cg2o build,
put G2O_USE_OPTIMIZATION_LIBRARY(umfpack_fast) into the program that uses the
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
#include "linear_solver_umfpack_fast.h"

namespace cg2o {

namespace {
template <int p, int l, bool blockorder>
std::unique_ptr<g2o::BlockSolverBase> AllocateSolverUmfpackFast() {
  G2O_DEBUG("Using LinearSolverUmfpackFast poseDim {} landMarkDim {} "
            "blockordering {}",
            p, l, blockorder);
  auto linearSolver = std::make_unique<LinearSolverUmfpackFast<
      typename g2o::BlockSolverPL<p, l>::PoseMatrixType>>();
  linearSolver->setBlockOrdering(blockorder);
  return std::make_unique<g2o::BlockSolverPL<p, l>>(std::move(linearSolver));
}

g2o::OptimizationAlgorithm *
createSolverUmfpackFast(const std::string &fullSolverName) {
  static const std::map<std::string,
                        std::function<std::unique_ptr<g2o::BlockSolverBase>()>>
      solver_factories{
          {"var", &AllocateSolverUmfpackFast<-1, -1, true>},
          {"fix3_2", &AllocateSolverUmfpackFast<3, 2, true>},
          {"fix6_3", &AllocateSolverUmfpackFast<6, 3, true>},
          {"fix7_3", &AllocateSolverUmfpackFast<7, 3, true>},
          {"var_scalar", &AllocateSolverUmfpackFast<-1, -1, false>},
          {"fix3_2_scalar", &AllocateSolverUmfpackFast<3, 2, false>},
          {"fix6_3_scalar", &AllocateSolverUmfpackFast<6, 3, false>},
          {"fix7_3_scalar", &AllocateSolverUmfpackFast<7, 3, false>},
      };
  // <method>_<blocksize>_umfpack_fast
  const std::string suffix = "_umfpack_fast";
  if (fullSolverName.size() <= 3 + suffix.size() ||
      fullSolverName.compare(fullSolverName.size() - suffix.size(),
                             suffix.size(), suffix) != 0)
    return nullptr;
  const std::string method = fullSolverName.substr(0, 2);
  const std::string blocks =
      fullSolverName.substr(3, fullSolverName.size() - 3 - suffix.size());
  auto solverf = solver_factories.find(blocks);
  if (solverf == solver_factories.end())
    return nullptr;
  if (method == "gn")
    return new g2o::OptimizationAlgorithmGaussNewton(solverf->second());
  if (method == "lm")
    return new g2o::OptimizationAlgorithmLevenberg(solverf->second());
  if (method == "dl")
    return new g2o::OptimizationAlgorithmDogleg(solverf->second());
  return nullptr;
}

class UmfpackFast_SolverCreator
    : public g2o::AbstractOptimizationAlgorithmCreator {
public:
  explicit UmfpackFast_SolverCreator(
      const g2o::OptimizationAlgorithmProperty &p)
      : g2o::AbstractOptimizationAlgorithmCreator(p) {}
  g2o::OptimizationAlgorithm *construct() override {
    return createSolverUmfpackFast(property().name);
  }
};
} // namespace

G2O_REGISTER_OPTIMIZATION_LIBRARY(umfpack_fast);

#define CG2O_REGISTER_UmfpackFast(method, blocks, desc, fixed, p, l)           \
  G2O_REGISTER_OPTIMIZATION_ALGORITHM(                                         \
      method##_##blocks##_umfpack_fast,                                        \
      new UmfpackFast_SolverCreator(g2o::OptimizationAlgorithmProperty(        \
          #method "_" #blocks "_umfpack_fast", desc, "UMFPACK", fixed, p, l)))

CG2O_REGISTER_UmfpackFast(gn, var,
                          "Gauss-Newton: UMFPACK LU with cached analysis "
                          "(SuiteSparse) (variable blocksize)",
                          false, Eigen::Dynamic, Eigen::Dynamic);
CG2O_REGISTER_UmfpackFast(gn, fix3_2,
                          "Gauss-Newton: UMFPACK LU with cached analysis "
                          "(SuiteSparse) (fixed blocksize)",
                          true, 3, 2);
CG2O_REGISTER_UmfpackFast(gn, fix6_3,
                          "Gauss-Newton: UMFPACK LU with cached analysis "
                          "(SuiteSparse) (fixed blocksize)",
                          true, 6, 3);
CG2O_REGISTER_UmfpackFast(gn, fix7_3,
                          "Gauss-Newton: UMFPACK LU with cached analysis "
                          "(SuiteSparse) (fixed blocksize)",
                          true, 7, 3);
CG2O_REGISTER_UmfpackFast(lm, var,
                          "Levenberg: UMFPACK LU with cached analysis "
                          "(SuiteSparse) (variable blocksize)",
                          false, Eigen::Dynamic, Eigen::Dynamic);
CG2O_REGISTER_UmfpackFast(lm, fix3_2,
                          "Levenberg: UMFPACK LU with cached analysis "
                          "(SuiteSparse) (fixed blocksize)",
                          true, 3, 2);
CG2O_REGISTER_UmfpackFast(lm, fix6_3,
                          "Levenberg: UMFPACK LU with cached analysis "
                          "(SuiteSparse) (fixed blocksize)",
                          true, 6, 3);
CG2O_REGISTER_UmfpackFast(lm, fix7_3,
                          "Levenberg: UMFPACK LU with cached analysis "
                          "(SuiteSparse) (fixed blocksize)",
                          true, 7, 3);
CG2O_REGISTER_UmfpackFast(dl, var,
                          "Dogleg: UMFPACK LU with cached analysis "
                          "(SuiteSparse) (variable blocksize)",
                          false, Eigen::Dynamic, Eigen::Dynamic);
CG2O_REGISTER_UmfpackFast(
    dl, fix3_2,
    "Dogleg: UMFPACK LU with cached analysis (SuiteSparse) (fixed blocksize)",
    true, 3, 2);
CG2O_REGISTER_UmfpackFast(
    dl, fix6_3,
    "Dogleg: UMFPACK LU with cached analysis (SuiteSparse) (fixed blocksize)",
    true, 6, 3);
CG2O_REGISTER_UmfpackFast(
    dl, fix7_3,
    "Dogleg: UMFPACK LU with cached analysis (SuiteSparse) (fixed blocksize)",
    true, 7, 3);

#undef CG2O_REGISTER_UmfpackFast

} // namespace cg2o
