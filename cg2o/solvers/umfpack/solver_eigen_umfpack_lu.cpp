#include "g2o/core/block_solver.h"
#include "g2o/core/optimization_algorithm_dogleg.h"
#include "g2o/core/optimization_algorithm_factory.h"
#include "g2o/core/optimization_algorithm_gauss_newton.h"
#include "g2o/core/optimization_algorithm_levenberg.h"
#include "g2o/core/sparse_optimizer.h"
#include "g2o/stuff/logger.h"
#include "linear_solver_eigen_umfpack_lu.h"  // Updated header

using namespace std;

namespace cg2o {

namespace {
template <int p, int l, bool blockorder>
std::unique_ptr<g2o::BlockSolverBase> AllocateSolverUmfPackLU() {
  G2O_DEBUG(
      "Using EigenUmfPackLU poseDim {} landMarkDim {} blockordering {}", p,
      l, blockorder);
  auto linearSolver = std::make_unique<
      LinearSolverEigenUmfPackLU<typename g2o::BlockSolverPL<p, l>::PoseMatrixType>>();
  linearSolver->setBlockOrdering(blockorder);
  return std::make_unique<g2o::BlockSolverPL<p, l>>(std::move(linearSolver));
}
}  // namespace

/**
 * Helper function for allocating UmfPackLU-based solvers
 */
static g2o::OptimizationAlgorithm* createSolverUmfPackLU(const std::string& fullSolverName) {
  static const std::map<std::string,
                        std::function<std::unique_ptr<g2o::BlockSolverBase>()>>
      solver_factories{
          {"var", &AllocateSolverUmfPackLU<-1, -1, true>},
          {"fix3_2", &AllocateSolverUmfPackLU<3, 2, true>},
          {"fix6_3", &AllocateSolverUmfPackLU<6, 3, true>},
          {"fix7_3", &AllocateSolverUmfPackLU<7, 3, true>},
          {"fix3_2_scalar", &AllocateSolverUmfPackLU<3, 2, false>},
          {"fix6_3_scalar", &AllocateSolverUmfPackLU<6, 3, false>},
          {"fix7_3_scalar", &AllocateSolverUmfPackLU<7, 3, false>},
      };

  string solverName = fullSolverName.substr(3);
  auto solverf = solver_factories.find(solverName);
  if (solverf == solver_factories.end()) return nullptr;

  string methodName = fullSolverName.substr(0, 2);

  if (methodName == "gn") {
    return new g2o::OptimizationAlgorithmGaussNewton(solverf->second());
  } else if (methodName == "lm") {
    return new g2o::OptimizationAlgorithmLevenberg(solverf->second());
  } else if (methodName == "dl") {
    return new g2o::OptimizationAlgorithmDogleg(solverf->second());
  }

  return nullptr;
}

class EigenUmfPackLU_SolverCreator : public g2o::AbstractOptimizationAlgorithmCreator {
 public:
  explicit EigenUmfPackLU_SolverCreator(const g2o::OptimizationAlgorithmProperty& p)
      : g2o::AbstractOptimizationAlgorithmCreator(p) {}
  virtual g2o::OptimizationAlgorithm* construct() {
    return createSolverUmfPackLU(property().name);
  }
};

// Register UmfPackLU-based solvers
G2O_REGISTER_OPTIMIZATION_LIBRARY(eigen_umfpack_lu);

G2O_REGISTER_OPTIMIZATION_ALGORITHM(gn_var_umfpack_lu, new EigenUmfPackLU_SolverCreator(
    g2o::OptimizationAlgorithmProperty("gn_var_umfpack_lu", "Gauss-Newton: UmfPackLU solver using SuiteSparse (variable blocksize)", "Eigen", false, Eigen::Dynamic, Eigen::Dynamic)));
G2O_REGISTER_OPTIMIZATION_ALGORITHM(gn_fix3_2_umfpack_lu, new EigenUmfPackLU_SolverCreator(
    g2o::OptimizationAlgorithmProperty("gn_fix3_2_umfpack_lu", "Gauss-Newton: UmfPackLU solver using SuiteSparse (fixed blocksize)", "Eigen", true, 3, 2)));
G2O_REGISTER_OPTIMIZATION_ALGORITHM(gn_fix6_3_umfpack_lu, new EigenUmfPackLU_SolverCreator(
    g2o::OptimizationAlgorithmProperty("gn_fix6_3_umfpack_lu", "Gauss-Newton: UmfPackLU solver using SuiteSparse (fixed blocksize)", "Eigen", true, 6, 3)));
G2O_REGISTER_OPTIMIZATION_ALGORITHM(gn_fix7_3_umfpack_lu, new EigenUmfPackLU_SolverCreator(
    g2o::OptimizationAlgorithmProperty("gn_fix7_3_umfpack_lu", "Gauss-Newton: UmfPackLU solver using SuiteSparse (fixed blocksize)", "Eigen", true, 7, 3)));

G2O_REGISTER_OPTIMIZATION_ALGORITHM(lm_var_umfpack_lu, new EigenUmfPackLU_SolverCreator(
    g2o::OptimizationAlgorithmProperty("lm_var_umfpack_lu", "Levenberg: UmfPackLU solver using SuiteSparse (variable blocksize)", "Eigen", false, Eigen::Dynamic, Eigen::Dynamic)));
G2O_REGISTER_OPTIMIZATION_ALGORITHM(lm_fix3_2_umfpack_lu, new EigenUmfPackLU_SolverCreator(
    g2o::OptimizationAlgorithmProperty("lm_fix3_2_umfpack_lu", "Levenberg: UmfPackLU solver using SuiteSparse (fixed blocksize)", "Eigen", true, 3, 2)));
G2O_REGISTER_OPTIMIZATION_ALGORITHM(lm_fix6_3_umfpack_lu, new EigenUmfPackLU_SolverCreator(
    g2o::OptimizationAlgorithmProperty("lm_fix6_3_umfpack_lu", "Levenberg: UmfPackLU solver using SuiteSparse (fixed blocksize)", "Eigen", true, 6, 3)));
G2O_REGISTER_OPTIMIZATION_ALGORITHM(lm_fix7_3_umfpack_lu, new EigenUmfPackLU_SolverCreator(
    g2o::OptimizationAlgorithmProperty("lm_fix7_3_umfpack_lu", "Levenberg: UmfPackLU solver using SuiteSparse (fixed blocksize)", "Eigen", true, 7, 3)));

G2O_REGISTER_OPTIMIZATION_ALGORITHM(dl_var_umfpack_lu, new EigenUmfPackLU_SolverCreator(
    g2o::OptimizationAlgorithmProperty("dl_var_umfpack_lu", "Dogleg: UmfPackLU solver using SuiteSparse (variable blocksize)", "Eigen", false, Eigen::Dynamic, Eigen::Dynamic)));

}  // namespace cg2o
