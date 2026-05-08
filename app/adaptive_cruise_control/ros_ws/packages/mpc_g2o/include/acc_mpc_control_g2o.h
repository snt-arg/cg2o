#ifndef ACC_MPC_CONTROL_G2O_H
#define ACC_MPC_CONTROL_G2O_H

#include "cg2o/solvers/eigen/linear_solver_eigen_svd.h"
#include "mpc_cg2o/mpc_formulation.h"
#include "mpc_cg2o/mpc_parameters.h"

#include "g2o/core/block_solver.h"
#include <g2o/core/optimization_algorithm_gauss_newton.h>

#include "cg2o/core/sparse_optimizer_ispd.h"
#include "cg2o/solvers/eigen/linear_solver_eigen_qr.h"
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
    switch (-1) {
    case 0:
      linearSolver = std::make_unique<
          cg2o::LinearSolverEigenSVD<g2o::BlockSolverX::PoseMatrixType>>();

      break;
    default:
      linearSolver = std::make_unique<cg2o::LinearSolverEigenUmfPackLU<
          g2o::BlockSolverX::PoseMatrixType>>();
      break;
    }

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

    std::cout << "Optimizer initialized with settings: " << std::endl;

    std::cout << "setKappaInitial: " << _optimizer->kappaInitial() << std::endl;
    std::cout << "setKappaFinal: " << _optimizer->kappaFinal() << std::endl;
    std::cout << "setKappaUpdateFactor: " << _optimizer->kappaUpdateFactor()
              << std::endl;
  }

  void initializeMPCGraph(int mpcHorizon) {
    _mpcParameters = std::make_shared<cg2o::mpc::MPCParameters>();
    _mpc =
        std::make_shared<cg2o::mpc::MPCFormulation<cg2o::SparseOptimizerISPD>>(
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
    //_mpc->setInitialGuess(false); // Set the initial guess
    _mpc->setFixedInitialGuess(10.0);
  }

  double min_of_vector(const std::vector<double> &data) {
    if (data.empty())
      return std::numeric_limits<double>::quiet_NaN();
    return *std::min_element(data.begin(), data.end());
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

    bool use_space_domain =
        _mpcParameters->isSpaceDomain() &&
        (min_of_vector(_mpcParameters->get_vp_prediction()) >
         _mpcParameters->get_domain_switch_threshold());

    std::cout << "min_of_vector(_mpcParameters->get_vp_prediction()): "
              << min_of_vector(_mpcParameters->get_vp_prediction())
              << std::endl;

    std::cout << "Domain switch threshold: "
              << _mpcParameters->get_domain_switch_threshold() << std::endl;

    _mpcParameters->set_use_space_domain(use_space_domain);
    if (use_space_domain) {
      std::cout
          << "*************************************************************"
             "*************************************************************"
             "******** SPACE DOMAIN "
          << std::endl;
    }

    if (_mpcParameters->isSpaceDomain()) {
      double v_p = prediction_data[0][0];
      v_p = v_p + 0.0;
      double delta_s =
          1.0 * _mpcParameters->get_v_h_0() * _mpcParameters->get_delta_t();
      _mpcParameters->set_delta_s(delta_s);
      double delta_s_min = 1.0;
      if (delta_s < delta_s_min) {
        _mpcParameters->set_delta_s(delta_s_min);
      }
    }
  }

  void printAll() {
    std::cout << "MPC parameters updated: " << std::endl;
    std::cout << "v_h_0: " << _mpcParameters->get_v_h_0() << std::endl;
    std::cout << "d_h_0: " << _mpcParameters->get_d_h_0() << std::endl;
    std::cout << "f_t_prev: " << _mpcParameters->get_f_t_prev() << std::endl;
    std::cout << "f_b_prev: " << _mpcParameters->get_f_b_prev() << std::endl;
    std::cout << "vp_prediction: "
              << Eigen::Map<Eigen::VectorXd>(
                     _mpcParameters->get_vp_prediction().data(),
                     _mpcParameters->get_vp_prediction().size())
                     .transpose()
              << std::endl;
    std::cout << "delta_s: " << _mpcParameters->get_delta_s() << std::endl;
    std::cout << "use_space_domain: " << _mpcParameters->useSpaceDomain()
              << std::endl;
    std::cout << " _param->isLinearInequalities() "
              << _mpcParameters->isLinearInequalities() << std::endl;
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
  std::string _solverType;
  std::string _algEqType;
  std::vector<double> _algIneqSettings;
  std::unique_ptr<g2o::OptimizationAlgorithm> _algorithm;
  std::shared_ptr<cg2o::SparseOptimizerISPD> _optimizer;
  std::shared_ptr<cg2o::mpc::MPCParameters> _mpcParameters;
  std::shared_ptr<cg2o::mpc::MPCFormulation<cg2o::SparseOptimizerISPD>> _mpc;
};

#endif // ACC_MPC_CONTROL_G2O_H
