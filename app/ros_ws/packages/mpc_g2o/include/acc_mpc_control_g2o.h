#ifndef ACC_MPC_CONTROL_G2O_H
#define ACC_MPC_CONTROL_G2O_H

#include "mpc_cg2o/mpc_formulation.h"
#include "mpc_cg2o/mpc_parameters.h"

#include "g2o/core/block_solver.h"
#include <g2o/core/optimization_algorithm_gauss_newton.h>

#include "cg2o/core/sparse_optimizer_ispd.h"
#include "cg2o/solvers/umfpack/linear_solver_eigen_umfpack_lu.h"

#include <iostream>
#include <memory>
#include <vector>
class ACCMPCControlG2O {
public:
  ACCMPCControlG2O(int numberOfIterations, std::vector<double> &AlgIneqSettings)
      : _numberOfIterations(numberOfIterations),        
        _algIneqSettings(AlgIneqSettings) {}

  void initializeOptimizer() {
    // Validate solver type

    std::unique_ptr<g2o::LinearSolver<g2o::BlockSolverX::PoseMatrixType>>
        linearSolver;
    linearSolver = std::make_unique<
        cg2o::LinearSolverEigenUmfPackLU<g2o::BlockSolverX::PoseMatrixType>>();

    auto blockSolver =
        std::make_unique<g2o::BlockSolverX>(std::move(linearSolver));

    _algorithm = std::make_unique<g2o::OptimizationAlgorithmGaussNewton>(
        std::move(blockSolver));

    _optimizer = std::make_shared<cg2o::SparseOptimizerISPD>();
    _optimizer->setVerbose(false);
    _optimizer->setAlgorithm(_algorithm.get());
    _optimizer->setKappaInitial(_algIneqSettings[0]);
    _optimizer->setKappaFinal(_algIneqSettings[1]);
    _optimizer->setKappaUpdateFactor(_algIneqSettings[2]);

   }

  void initializeMPCGraph(int mpcHorizon) {
    _mpcParameters = std::make_shared<cg2o::mpc::MPCParameters>();
    _mpc = std::make_shared<cg2o::mpc::MPCFormulation<cg2o::SparseOptimizerISPD>>(
        mpcHorizon, _optimizer, _mpcParameters);
    _mpc->setupMPC();
    _optimizer->initializeOptimization();
  }

  void updateMPCParameters(double v_h, double d_h, double f_t_prev,
                           double f_b_prev) {
    _mpcParameters->set_v_h_0(v_h);
    _mpcParameters->set_d_h_0(d_h);
    _mpcParameters->set_f_t_prev(f_t_prev);
    _mpcParameters->set_f_b_prev(f_b_prev);

    _mpc->setInitialGuess(false);
    _mpc->setFixedInitialGuess(10.0);
  }

  void
  updateMPCPrediction(const std::vector<std::vector<double>> &prediction_data) {
    if (prediction_data.size() < 4) {
      throw std::invalid_argument(
          "Prediction data must have at least 4 vectors.");
    }
    _mpcParameters->set_vp_prediction(prediction_data[0]);
    _mpcParameters->set_alpha_prediction(prediction_data[1]);
    _mpcParameters->set_v_min_prediction(prediction_data[2]);
    _mpcParameters->set_v_max_prediction(prediction_data[3]);
  }

  void performOptimization() {
    if (!_optimizer) {
      throw std::runtime_error(
          "Optimizer not initialized. Call initializeOptimizer first.");
    }

    auto results = _mpc->getResults();
    std::cout << "Final Resutls:" << std::endl;
    std::cout << Eigen::Map<Eigen::VectorXd>(results.data(), results.size())
                     .transpose()
              << std::endl;
     
    
    _optimizer->optimize(_numberOfIterations);

    results = _mpc->getResults();
    double u_input = _mpc->getForceInput(0);
    std::cout << "Optimization done.\nResults:\n"
              << Eigen::Map<Eigen::VectorXd>(results.data(), results.size())
                     .transpose()
              << std::endl;
    std::cout << "u_input: " << u_input << std::endl;
  }

  void resetGraph() {
    if (_optimizer) {
      _optimizer->clear();
    }
  }

  void setNumberOfIterations(int numberOfIterations) {
    _numberOfIterations = numberOfIterations;
  }

  int getNumberOfIterations() const { return _numberOfIterations; }

  std::shared_ptr<cg2o::SparseOptimizerISPD> optimizer() const {
    return _optimizer;
  }

  std::shared_ptr<cg2o::mpc::MPCFormulation<cg2o::SparseOptimizerISPD>>
  mpc() const {
    return _mpc;
  }

private:
  int _numberOfIterations;
  int _linearSolverType;
  std::string _solverType;
  std::string _algEqType;
   std::vector<double> _algIneqSettings;
   std::unique_ptr<g2o::OptimizationAlgorithm> _algorithm;
  std::shared_ptr<cg2o::SparseOptimizerISPD> _optimizer;
  std::shared_ptr<cg2o::mpc::MPCParameters> _mpcParameters;
  std::shared_ptr<cg2o::mpc::MPCFormulation<cg2o::SparseOptimizerISPD>> _mpc;
};

#endif // ACC_MPC_CONTROL_G2O_H
