#include "sparse_optimizer_ispd.h"
#include "g2o/core/optimization_algorithm_with_hessian.h"
#include "g2o/core/sparse_optimizer.h"
#include "g2o/stuff/logger.h"
#include "g2o/stuff/timeutil.h"
#include "optimization_algorithm.h"
#include "solver.h"
#include <Eigen/Core>
#include <cmath>
#include <iostream>
#include <memory>
#include <vector>

namespace cg2o {
using namespace std;

// Default constructor
SparseOptimizerISPD::SparseOptimizerISPD() = default;

// Default destructor
SparseOptimizerISPD::~SparseOptimizerISPD() = default;

// Setter solver paramters
void SparseOptimizerISPD::setStepSizeStrategy(int strategy) {
  _step_size_strategy = strategy;
}

void SparseOptimizerISPD::setIneqBacktrackingStepMin(double value) {
  _ineq_backtracking_step_min = value;
}

void SparseOptimizerISPD::setAuxBacktrackingStepMin(double value) {
  _aux_backtracking_step_min = value;
}

void SparseOptimizerISPD::setInitKappaStrategy(int strategy) {
  _init_kappa_strategy = strategy;
}

void SparseOptimizerISPD::setKappaInitial(double value) {
  _kappa_initial = value;
}

void SparseOptimizerISPD::setUpdateKappaStrategy(int strategy) {
  _update_kappa_strategy = strategy;
}

void SparseOptimizerISPD::setKappaFinal(double value) { _kappa_final = value; }

void SparseOptimizerISPD::setKappaUpdateFactor(double value) {
  _kappa_update_factor = value;
}

void SparseOptimizerISPD::setTau(double value) { _tau = value; }

void SparseOptimizerISPD::setLimitKappaFinal(bool value) {
  _limit_kappa_final = value;
}

void SparseOptimizerISPD::setInitSlackStrategy(int strategy) {
  _init_slack_strategy = strategy;
}

void SparseOptimizerISPD::setSlackVariableInitialIneq(double value) {
  _slack_variable_initial_ineq = value;
}

void SparseOptimizerISPD::setInitLagrangeStrategy(int strategy) {
  _init_lagrange_strategy = strategy;
}

void SparseOptimizerISPD::setLagrangeMultiplierInitialIneq(double value) {
  _lagrange_multiplier_initial_ineq = value;
}

void SparseOptimizerISPD::setIneqPredictionStrategy(int strategy) {
  _ineq_prediction_strategy = strategy;
}

void SparseOptimizerISPD::setIneqPredictionStepStrategy(int strategy) {
  _ineq_prediction_step_strategy = strategy;
}

void SparseOptimizerISPD::setIneqPredictionParam(double value) {
  _ineq_prediction_param = value;
}

void SparseOptimizerISPD::setKeepAuxPositiveStrategy(int strategy) {
  _keep_aux_positive_strategy = strategy;
}

void SparseOptimizerISPD::setAuxScalingFactor(double value) {
  _aux_scaling_factor = value;
}

void SparseOptimizerISPD::setAuxCorrectionValue(double value) {
  _aux_correction_value = value;
}

int SparseOptimizerISPD::stepSizeStrategy() const {
  return _step_size_strategy;
}

double SparseOptimizerISPD::ineqBacktrackingStepMin() const {
  return _ineq_backtracking_step_min;
}

double SparseOptimizerISPD::auxBacktrackingStepMin() const {
  return _aux_backtracking_step_min;
}

int SparseOptimizerISPD::initKappaStrategy() const {
  return _init_kappa_strategy;
}

double SparseOptimizerISPD::kappaInitial() const { return _kappa_initial; }

int SparseOptimizerISPD::updateKappaStrategy() const {
  return _update_kappa_strategy;
}

double SparseOptimizerISPD::kappaUpdateFactor() const {
  return _kappa_update_factor;
}

double SparseOptimizerISPD::kappaFinal() const { return _kappa_final; }

double SparseOptimizerISPD::tau() const { return _tau; }

bool SparseOptimizerISPD::limitKappaFinal() const { return _limit_kappa_final; }

int SparseOptimizerISPD::initSlackStrategy() const {
  return _init_slack_strategy;
}

double SparseOptimizerISPD::slackVariableInitialIneq() const {
  return _slack_variable_initial_ineq;
}

int SparseOptimizerISPD::initLagrangeStrategy() const {
  return _init_lagrange_strategy;
}

double SparseOptimizerISPD::lagrangeMultiplierInitialIneq() const {
  return _lagrange_multiplier_initial_ineq;
}

int SparseOptimizerISPD::ineqPredictionStrategy() const {
  return _ineq_prediction_strategy;
}

int SparseOptimizerISPD::ineqPredictionStepStrategy() const {
  return _ineq_prediction_step_strategy;
}

double SparseOptimizerISPD::ineqPredictionParam() const {
  return _ineq_prediction_param;
}

int SparseOptimizerISPD::keepAuxPositiveStrategy() const {
  return _keep_aux_positive_strategy;
}

double SparseOptimizerISPD::auxScalingFactor() const {
  return _aux_scaling_factor;
}

double SparseOptimizerISPD::auxCorrectionValue() const {
  return _aux_correction_value;
}

void SparseOptimizerISPD::setAlphaBacktracking(
    std::vector<double> alphaBacktracking) {
  _alphaBacktracking = alphaBacktracking;
}

void SparseOptimizerISPD::resetLagrangeMultiplierEq() {
  for (auto &vertex : _vEqLagrangeMultipliers) {
    vertex->setToOrigin();
  }
}

void SparseOptimizerISPD::executeEdgeProcessing(void *edgePtr, int controller) {
  auto it = edgeProcessingFunctionMap.find(edgePtr);
  if (it != edgeProcessingFunctionMap.end()) {
    it->second(controller); // Calls updateMultipliers on the correct edge
  } else {
    std::cerr << "[Error] Edge not found in function map in Update function"
              << std::endl;
  }
}

void SparseOptimizerISPD::computeKappa(bool input) {
  // limit t to t_final
  if (_limit_kappa_final && (_t >= _kappa_final)) {
    return;
  }

  if (_init_kappa_strategy == 3) {
    if (input) {
      _t = _kappa_update_factor * _t;
    }
    return;
  }

  double eta = 0.0;
  int edges_counter = 0;
  for (auto &edge : _activeEdgesIneq) {
    auto it = getDualityGapFunctionMap.find(edge);
    if (it != getDualityGapFunctionMap.end()) {
      eta += it->second(); // Calls the function to compute the duality gap
      edges_counter = edges_counter + edge->dimension();
    } else {
      std::cerr << "[Error] Edge not found in function map" << std::endl;
    }
  }

  // ... compute eta, edges_counter

  double t_candidate = edges_counter * _kappa_update_factor / eta;

  switch (_update_kappa_strategy) {
  case 0:

    t_candidate = std::max(t_candidate,
                           (_tau + std::exp(-_t / (.2 * _kappa_final))) * _t);
    break;
  case 1:
    t_candidate = std::max(t_candidate, _kappa_initial);
    break;
  case 2:
    // keep t_candidate
    break;
  case 3:
    t_candidate = _t;
    if (input) {
      t_candidate = std::max(t_candidate, _kappa_initial);
      t_candidate = _t * _kappa_update_factor;
    }
    break;
  default:
    throw std::runtime_error("[Error] Unknown t update strategy.");
  }

  _t = t_candidate;
}

double SparseOptimizerISPD::constraintsBacktracking(const double *update) {
  int w = 0;
  // Initialize variables
  size_t c = 0;
  std::vector<bool> _inside_interior_flag;
  // Loop over all inequality edges
  for (auto &edge : _activeEdgesIneq) {
    w = w + 1;
    _inside_interior_flag.resize(edge->dimension());
    edge->computeError();
    // check if all is feasible
    for (int i = 0; i < edge->dimension(); ++i) {
      if (edge->errorData()[i] > 0) {
        _inside_interior_flag[i] = false;
      } else {
        _inside_interior_flag[i] = true;
      }
    }

    while (c < _alphaBacktracking.size() - 1 &&
           _alphaBacktracking[c] > _ineq_backtracking_step_min) {
      bool backtracking_done = true;
      // 1- Scale the update and apply it to the vertices
      updateEdgeVertices(edge, update, _alphaBacktracking[c]);

      // 2- Check if the inequality constraint is satisfied
      edge->computeError();
      auto error = edge->errorData();

      // 3- Reverse the update and apply it to the vertices
      updateEdgeVertices(edge, update, -_alphaBacktracking[c]);

      // 4- Check if the inequality constraint is satisfied
      for (int i = 0; i < edge->dimension(); ++i) {
        if (_inside_interior_flag[i] && error[i] > 0) {
          backtracking_done = false;
          break;
        }
      }

      // 5- Break if the inequality constraint is satisfied
      if (backtracking_done)
        break;

      // 6- Increment the step size index if the inequality constraint is not
      // satisfied
      ++c;
    }
  }
  return std::max(_alphaBacktracking[c], _ineq_backtracking_step_min);
}

void SparseOptimizerISPD::update(const double *update) {
  this->_update_step = const_cast<double *>(update);
  double ineq_backtracking_step = 1.0;
  double aux_backtracking_step = 1.0;
  // 1.a backtraking for inequality constraints
  if (_ineq_backtracking_step_min < 1) {
    ineq_backtracking_step = constraintsBacktracking(update);
  }

  // 1.b backtraking for auxiliary variables
  if (_aux_backtracking_step_min < 1) {
    _step_size = ineq_backtracking_step; // for the prediction step
    _step_size_index = 0;

    for (auto &edge : _activeEdgesIneq) {
      // controller =1: compute the error prediction
      executeEdgeProcessing(edge, 1);
    }

    for (auto &edge : _activeEdgesIneq) {
      // controller =2: backtracking the lagrange multiplier and slack variable
      executeEdgeProcessing(edge, 2);
    }
    aux_backtracking_step = std::max(_alphaBacktracking[_step_size_index],
                                     _aux_backtracking_step_min);
  }

  // 2.a - update the lagrange multiplier and slack variable using
  // aux_backtracking_step
  _step_size =
      aux_backtracking_step; // needed to compute the ineq prediction step if
                             // the strategy is based on the step size
  for (auto &edge : _activeEdgesIneq) {
    // controller =3: update the lagrange multiplier and slack variable
    executeEdgeProcessing(edge, 3);
  }
  switch (_step_size_strategy) {
  case 0: {
    _step_size = std::min(ineq_backtracking_step, aux_backtracking_step);
    break;
  }
  case 1: {
    _step_size = (ineq_backtracking_step + aux_backtracking_step) / 2.0;
    break;
  }
  case 2: {
    _step_size = ineq_backtracking_step;
    break;
  }
  case 3: {
    _step_size = aux_backtracking_step;
    break;
  }
  case 4: {
    _step_size = 1.0;
    break;
  }
  default:
    throw std::runtime_error("[Error] Unknown step size strategy");
  }
  // 2.b- update the decision variables
  g2o::OptimizationAlgorithmWithHessian *algorithm =
      static_cast<g2o::OptimizationAlgorithmWithHessian *>(_algorithm);
  size_t xSize = algorithm->solver().vectorSize();
  Eigen::Map<const Eigen::VectorXd> updateVec(update, xSize);
  Eigen::VectorXd scaledUpdateVec;
  const double *requiredUpdate = nullptr;
  if (_step_size == 1) { // No scaling required
    requiredUpdate = update;
  } else {
    scaledUpdateVec = _step_size * updateVec;
    requiredUpdate = scaledUpdateVec.data();
  }

  // std::cout << "....... ineq_step: " << ineq_backtracking_step
  //          << ", aux_step: " << aux_backtracking_step << ", applied_step "
  //<< _step_size << std::endl;

  SparseOptimizer::update(requiredUpdate);
}

int SparseOptimizerISPD::optimize(int iterations, bool online) {
  using namespace g2o;
  if (_ivMap.size() == 0) {
    G2O_WARN("0 vertices to optimize, maybe forgot to call "
             "initializeOptimization()");
    return -1;
  }

  int cjIterations = 0;
  double cumTime = 0;
  bool ok = true;

  bool stop = false;
  resetLagrangeMultiplierEq();

  if (_activeEdgesIneq.empty()) {
    _t = _kappa_final;
  } else {
    _t = _kappa_initial; // compute the initial value of _t or you can to
                         // compute t using
  }

  for (auto &edge : _activeEdgesIneq) {
    // initialize the lagrange multipliers and slack variables
    executeEdgeProcessing(edge, 0);
  }

  if (_init_kappa_strategy == 1) {
    computeKappa(false);
  }

  ok = _algorithm->init(online);
  if (!ok) {
    G2O_ERROR("Error while initializing");
    return -1;
  }

  _batchStatistics.clear();
  if (_computeBatchStatistics)
    _batchStatistics.resize(iterations);

  OptimizationAlgorithm::SolverResult result = OptimizationAlgorithm::OK;

  int i = 0;
  bool kappa_flag = false;

  while (!stop) {

    preIteration(cjIterations);
    if (_computeBatchStatistics) {
      G2OBatchStatistics &cstat = _batchStatistics[cjIterations];
      G2OBatchStatistics::setGlobalStats(&cstat);
      cstat.iteration = cjIterations;
      cstat.numEdges = _activeEdges.size();
      cstat.numVertices = _activeVertices.size();
    }

    double ts = get_monotonic_time();
    result = _algorithm->solve(cjIterations, online);
    ok = (result == OptimizationAlgorithm::OK);

    bool errorComputed = false;
    if (_computeBatchStatistics) {
      computeActiveErrors();
      errorComputed = true;
      _batchStatistics[cjIterations].chi2 = activeRobustChi2();
      _batchStatistics[cjIterations].timeIteration = get_monotonic_time() - ts;
    }

    if (verbose()) {
      double dts = get_monotonic_time() - ts;
      cumTime += dts;
      if (!errorComputed)
        computeActiveErrors();
      cerr << "iteration= " << cjIterations
           << "\t chi2= " << FIXED(activeRobustChi2()) << "\t time= " << dts
           << "\t cumTime= " << cumTime << "\t edges= " << _activeEdges.size();
      _algorithm->printVerbose(cerr);
      cerr << endl;
    }
    // std::cout << ", " << cjIterations << ": " << _t << " "
    //       << FIXED(activeRobustChi2());
    ++cjIterations;
    postIteration(cjIterations);
    // std::cout << "*************** _t: " << _t << ", cjIterations: *"
    //      << cjIterations << std::endl;

    stop = cjIterations >= iterations || terminate() || !ok;

    // termination criteria
    if (_update_kappa_strategy == 3) {
      i++;
      bool is_loop_done = verifyConvergence(-1 * 10) || stop ||
                          (i >= _num_inner_iterations_max);

      kappa_flag = is_loop_done;

      if (is_loop_done) {
        i = 0;
      }
    }

    if (cjIterations >= iterations ||
        terminate()) { // check if the maximum number of iterations is reached
      break;
    }

    if (result == OptimizationAlgorithm::Fail) {
      return 0;
    }

    bool converged =
        verifyConvergence(-1) && verifyEqFeasibility(_activeEdgesEq, -1.0) &&
        verifyIneqFeasibility(_activeEdgesIneq, -1.0) && _t >= _kappa_final;

    stop = stop || converged;

    if (!stop) {
      computeKappa(kappa_flag);
    }
  }
  // std::cout <<
  // "***********************************************************************************************************************************************************************"
  //          << cjIterations << std::endl;
  std::cout << std::endl;
  return cjIterations;
}

} // namespace cg2o
