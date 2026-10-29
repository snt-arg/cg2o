#include "sparse_optimizer_ispd.h"
#include "g2o/core/optimization_algorithm_with_hessian.h"
#include "g2o/core/sparse_optimizer.h"
#include "g2o/stuff/logger.h"
#include "g2o/stuff/timeutil.h"
#include "optimization_algorithm.h"
#include "solver.h"
#include "vertex_state_guard.h"
#include "warm_start_math.h"
#include <Eigen/Core>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <vector>

namespace cg2o {
using namespace std;

// Default constructor
SparseOptimizerISPD::SparseOptimizerISPD() = default;

// Default destructor
SparseOptimizerISPD::~SparseOptimizerISPD() = default;

bool SparseOptimizerISPD::removeEdge(g2o::HyperGraph::Edge *e) {
  if (!edges().count(e))
    return false;
  auto *oe = static_cast<g2o::OptimizableGraph::Edge *>(e);
  _start_policy.forgetEdge(oe);
  edgeProcessingFunctionMap.erase(oe);
  getDualityGapFunctionMap.erase(oe);
  _ineq_blocks.erase(oe);
  return SparseOptimizerForConstraints<SparseOptimizerISPD>::removeEdge(e);
}

bool SparseOptimizerISPD::removeVertex(g2o::HyperGraph::Vertex *v,
                                       bool detach) {
  if (!v || vertex(v->id()) != v)
    return false;
  _start_policy.forgetVertex(static_cast<g2o::OptimizableGraph::Vertex *>(v));
  return SparseOptimizerForConstraints<SparseOptimizerISPD>::removeVertex(
      v, detach);
}

void SparseOptimizerISPD::clear() {
  _start_policy.clear();
  edgeProcessingFunctionMap.clear();
  getDualityGapFunctionMap.clear();
  _ineq_blocks.clear();
  _last_converged = false;
  SparseOptimizerForConstraints<SparseOptimizerISPD>::clear();
}

void SparseOptimizerISPD::validateSettings() const {
  auto positive = [](double x) { return warm_start::positiveFinite(x); };
  if (!positive(_kappa_initial) || !positive(_kappa_final) ||
      !positive(1.0 / _kappa_initial) || !positive(1.0 / _kappa_final) ||
      (_limit_kappa_final && _kappa_final < _kappa_initial) ||
      !positive(_kappa_update_factor) ||
      !positive(_slack_variable_initial_ineq) ||
      !positive(_lagrange_multiplier_initial_ineq))
    throw std::invalid_argument(
        "ISPD needs finite positive kappa and auxiliary initialization; "
        "kappa_final >= kappa_initial");
  // an unknown init strategy would silently do nothing (lambda) or act as 0
  // (kappa); the others would throw in the middle of a solve
  auto inRange = [](int value, int last) {
    return value >= 0 && value <= last;
  };
  if (!inRange(_step_size_strategy, 4) || !inRange(_init_kappa_strategy, 1) ||
      !inRange(_update_kappa_strategy, 3) ||
      !inRange(_init_slack_strategy, 2) ||
      !inRange(_init_lagrange_strategy, 1) ||
      !inRange(_ineq_prediction_strategy, 1) ||
      !inRange(_ineq_prediction_step_strategy, 3) ||
      !inRange(_keep_aux_positive_strategy, 2))
    throw std::invalid_argument(
        "ISPD strategy out of range: step size 0-4, init kappa 0-1, update "
        "kappa 0-3, init slack 0-2, init lagrange 0-1, prediction 0-1, "
        "prediction step 0-3, keep aux positive 0-2");
  _start_policy.validate();
}

void SparseOptimizerISPD::resetLagrangeMultiplierEq() {
  for (auto &vertex : _vEqLagrangeMultipliers) {
    vertex->setToOrigin();
  }
}

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

// Getter solver parameters
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

const SparseOptimizerISPD::PointStats &
SparseOptimizerISPD::lastStartStats() const {
  return _start_stats;
}

const SparseOptimizerISPD::PointStats &
SparseOptimizerISPD::lastFirstIterationStats() const {
  return _first_stats;
}

void SparseOptimizerISPD::setRecordIterationStats(bool record) {
  _record_stats = record;
}

const std::vector<SparseOptimizerISPD::PointStats> &
SparseOptimizerISPD::iterationStats() const {
  return _iteration_stats;
}

SparseOptimizerISPD::PointStats SparseOptimizerISPD::pointStats(bool withCost) {
  PointStats st;
  double sum = 0.0;
  double smallest = std::numeric_limits<double>::infinity();
  double largest = -std::numeric_limits<double>::infinity();
  int rows = 0;
  forEachActiveIneqBlock([&](IneqBlock &block) {
    for (int i = 0; i < block.dimension; ++i) {
      const double p = block.slack[i] * block.multiplier[i];
      sum += p;
      smallest = std::min(smallest, p);
      largest = std::max(largest, p);
      ++rows;
    }
  });
  st.rows = rows;
  st.kappa = _kappa;
  if (rows > 0) {
    st.mu = sum / rows;
    if (st.mu > 0.0) {
      st.minRatio = smallest / st.mu;
      st.maxRatio = largest / st.mu;
      st.centrality =
          std::max(std::abs(st.minRatio - 1.0), std::abs(st.maxRatio - 1.0));
    }
  }
  if (withCost) {
    // Products and objective must describe the SAME iterate. Constraint
    // chi2() is zero; the robust objective also supports robust cost edges.
    computeActiveErrors();
    st.cost = 0.5 * activeRobustChi2();
    st.relativeGap = sum / (1.0 + std::abs(st.cost));
  }
  return st;
}

bool SparseOptimizerISPD::lastRestored() const {
  return _start_policy.lastRestored();
}

bool SparseOptimizerISPD::lastRetried() const {
  return _start_policy.lastRetried();
}
int SparseOptimizerISPD::lastResetRows() const {
  return _start_policy.lastResetRows();
}

void SparseOptimizerISPD::clearWarmStartShiftSequences() {
  _start_policy.clearShiftSequences();
}

void SparseOptimizerISPD::resetWarmStart() { _start_policy.invalidate(); }

bool SparseOptimizerISPD::warmStartAvailable() const {
  return _start_policy.available(_last_converged);
}

int SparseOptimizerISPD::lastWarmStartUsed() const {
  return _start_policy.lastSource();
}

double SparseOptimizerISPD::lastKappa() const { return _kappa; }

double SparseOptimizerISPD::lastStartKappa() const { return _last_start_kappa; }

bool SparseOptimizerISPD::lastConverged() const { return _last_converged; }

bool SparseOptimizerISPD::auxiliariesValid() {
  bool valid = true;
  forEachActiveIneqBlock([&](IneqBlock &block) {
    for (int i = 0; i < block.dimension; ++i)
      valid = valid && warm_start::validPair(block.slack[i],
                                             block.multiplier[i], _kappa);
  });
  return valid;
}

double SparseOptimizerISPD::dualityGapKappa() {
  double eta = 0.0;
  int m = 0;
  for (auto &edge : _activeEdgesIneq) {
    auto it = getDualityGapFunctionMap.find(edge);
    if (it != getDualityGapFunctionMap.end()) {
      eta += it->second();
      m += edge->dimension();
    }
  }
  if (m == 0 || !(eta > 0.0)) {
    return std::numeric_limits<double>::infinity();
  }
  return static_cast<double>(m) / eta;
}

double SparseOptimizerISPD::constraintViolation() {
  double violation = 0.0;
  for (auto &edge : _activeEdgesEq) {
    edge->computeError();
    const double *error = edge->errorData();
    const int dim = edge->dimension() / 2; // [h; nu]
    for (int i = 0; i < dim; ++i) {
      if (!std::isfinite(error[i])) {
        return std::numeric_limits<double>::infinity();
      }
      violation = std::max(violation, std::abs(error[i]));
    }
  }
  for (auto &edge : _activeEdgesIneq) {
    edge->computeError();
    const double *error = edge->errorData();
    for (int i = 0; i < edge->dimension(); ++i) {
      if (!std::isfinite(error[i])) {
        return std::numeric_limits<double>::infinity();
      }
      violation = std::max(violation, error[i]);
    }
  }
  return violation;
}

void SparseOptimizerISPD::coldStart(bool forceEqualityReset) {
  // keepEqualityMultipliers (setWarmStartLagrangeMultiplierFlag) keeps the
  // equality multipliers of the previous solve
  if (forceEqualityReset || !_start_policy.options().keepEqualityMultipliers) {
    resetLagrangeMultiplierEq();
  }

  if (_activeEdgesIneq.empty()) {
    _kappa = _kappa_final;
  } else {
    _kappa = _kappa_initial;
  }

  for (auto &edge : _activeEdgesIneq) {
    // initialize the lagrange multipliers and slack variables
    executeEdgeProcessing(edge, IneqStateOperation::ColdInitialize);
  }

  if (_init_kappa_strategy == 1) {
    computeKappa(false);
  }
}

// solver functions

void SparseOptimizerISPD::executeEdgeProcessing(void *edgePtr,
                                                IneqStateOperation operation) {
  auto it = edgeProcessingFunctionMap.find(edgePtr);
  if (it != edgeProcessingFunctionMap.end()) {
    it->second(operation); // Calls updateMultipliers on the correct edge
  } else {
    std::cerr << "[Error] Edge not found in function map in Update function"
              << std::endl;
  }
}

void SparseOptimizerISPD::computeKappa(bool input) {
  // limit t to t_final
  if (_limit_kappa_final && (_kappa >= _kappa_final)) {
    return;
  }

  // fixed-factor update (strategy 3); this tested _init_kappa_strategy
  // before, which froze kappa when init_kappa_strategy == 3
  if (_update_kappa_strategy == 3) {
    if (input) {
      _kappa = _kappa_update_factor * _kappa;
      if (_limit_kappa_final)
        _kappa = std::min(_kappa, _kappa_final);
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

  if (edges_counter == 0) {
    _kappa = _kappa_final;
    return;
  }
  if (!(eta > 0.0) || !std::isfinite(eta))
    return;
  double t_candidate = edges_counter * _kappa_update_factor / eta;

  switch (_update_kappa_strategy) {
  case 0:

    t_candidate = std::max(
        t_candidate, (_tau + std::exp(-_kappa / (.2 * _kappa_final))) * _kappa);
    break;
  case 1:
    t_candidate = std::max(t_candidate, _kappa_initial);
    break;
  case 2:
    // keep t_candidate
    break;
  case 3:
    t_candidate = _kappa;
    if (input) {
      t_candidate = std::max(t_candidate, _kappa_initial);
      t_candidate = _kappa * _kappa_update_factor;
    }
    break;
  default:
    throw std::runtime_error("[Error] Unknown t update strategy.");
  }

  if (_limit_kappa_final)
    t_candidate = std::min(t_candidate, _kappa_final);
  if (warm_start::positiveFinite(t_candidate))
    _kappa = t_candidate;
}

double SparseOptimizerISPD::constraintsBacktracking(const double *update) {
  // Initialize variables
  size_t c = 0;
  std::vector<bool> _inside_interior_flag;
  // Loop over all inequality edges
  for (auto &edge : _activeEdgesIneq) {
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
      internal::VertexStateGuard trial(edge->vertices());
      updateEdgeVertices(edge, update, _alphaBacktracking[c]);

      // 2- Check if the inequality constraint is satisfied
      edge->computeError();
      const double *error = edge->errorData();
      // Inspect the trial before restoring; no temporary error-vector copy.
      for (int i = 0; i < edge->dimension(); ++i) {
        if (!std::isfinite(error[i]) ||
            (_inside_interior_flag[i] && error[i] > 0)) {
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
    edge->computeError(); // refresh once, after all trials for this edge
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

  // 1.b predicted inequality values and the multiplier / slack updates,
  // with backtracking for the auxiliary variables. The updates are needed
  // even when the backtracking is switched off (_aux_backtracking_step_min
  // >= 1): before, they were then never computed and the update applied
  // stale values.
  _step_size = ineq_backtracking_step; // for the prediction step
  _step_size_index = 0;

  for (auto &edge : _activeEdgesIneq) {
    // compute the error prediction
    executeEdgeProcessing(edge, IneqStateOperation::PredictConstraint);
  }

  for (auto &edge : _activeEdgesIneq) {
    // updates, and backtracking the lagrange multiplier and
    // slack variable (no backtracking when _aux_backtracking_step_min >= 1)
    executeEdgeProcessing(edge, IneqStateOperation::ComputeAuxiliaryStep);
  }
  if (_aux_backtracking_step_min < 1) {
    aux_backtracking_step = std::max(_alphaBacktracking[_step_size_index],
                                     _aux_backtracking_step_min);
  }

  // 2.a - update the lagrange multiplier and slack variable using
  // aux_backtracking_step
  _step_size =
      aux_backtracking_step; // needed to compute the ineq prediction step if
                             // the strategy is based on the step size
  _scaled_count = 0;
  for (auto &edge : _activeEdgesIneq) {
    // update the lagrange multiplier and slack variable
    executeEdgeProcessing(edge, IneqStateOperation::ApplyAuxiliaryStep);
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
  _last_ineq_step = ineq_backtracking_step;
  _last_aux_step = aux_backtracking_step;
  _last_applied_step = _step_size;
  applyScaledUpdate(update, _step_size);
}

int SparseOptimizerISPD::optimize(int iterations, bool online) {
  using namespace g2o;
  BatchStatisticsScope statisticsScope;
  _last_iteration_count = 0;
  if (iterations < 0)
    throw std::invalid_argument("iteration budget must be >= 0");
  if (iterations == 0 || terminate())
    return 0;
  if (!fixedConstraintsFeasible()) {
    G2O_WARN("constraints on fixed vertices are violated, nothing to "
             "optimize");
    _last_converged = false;
    _start_policy.skipSolve();
    return 0;
  }
  validateSettings();
  if (_ivMap.size() == 0) {
    G2O_WARN("0 vertices to optimize, maybe forgot to call "
             "initializeOptimization()");
    return -1;
  }

  // starting point: cold (nu, lambda, s, kappa re-initialised) or warm
  // (cached point of the previous solve, optionally shifted), see
  // ISPDStartOptions::source; uses _last_converged of the previous solve
  const StartResult start = _start_policy.prepare(*this);
  _last_converged = false;
  _last_start_kappa = _kappa;
  _iteration_stats.clear();
  _start_stats = PointStats();
  _first_stats = PointStats();
  if (_start_policy.options().source != ISPDStartSource::Cold ||
      _record_stats) {
    // cost needs the errors at the starting point (one extra evaluation)
    _start_stats = pointStats(_record_stats);
    _start_stats.iteration = -1;
  }
  _start_stats.violation = start.startViolation;

  SolveOutcome outcome = runIterations(iterations, online);
  _last_iteration_count = std::max(outcome.iterations, 0);
  int done = outcome.iterations;
  FinishAction action = _start_policy.finish(*this, outcome);
  if (action == FinishAction::RetryCold) {
    _start_policy.prepareColdRetry(*this);
    outcome = runIterations(iterations, online);
    _last_iteration_count += std::max(outcome.iterations, 0);
    if (outcome.iterations >= 0)
      done += outcome.iterations;
    action = _start_policy.finish(*this, outcome);
  }
  switch (action) {
  case FinishAction::Error:
    return -1;
  case FinishAction::Reject:
    return outcome.solverFailed ? 0 : done; // 0 reports a solver failure
  default:
    return done;
  }
}

SolveOutcome SparseOptimizerISPD::runIterations(int iterations, bool online) {
  using namespace g2o;
  int cjIterations = 0;
  double cumTime = 0;
  bool ok = true;

  bool stop = false;

  ok = _algorithm->init(online);
  if (!ok) {
    G2O_ERROR("Error while initializing");
    return {-1, false};
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

    // statistics of the new iterate (first iteration: step sizes)
    const bool retry = _start_policy.lastRetried();
    const bool middle = ok && _start_policy.observesIterations();
    const bool needStats =
        _record_stats || middle ||
        (cjIterations == 0 &&
         _start_policy.options().source != ISPDStartSource::Cold);
    PointStats st;
    if (needStats) {
      st = pointStats(_record_stats ||
                      (middle && _start_policy.needsIterationCost()));
      st.iteration = cjIterations;
      st.attempt = retry ? 1 : 0;
      st.ineqStep = _last_ineq_step;
      st.auxStep = _last_aux_step;
      st.appliedStep = _last_applied_step;
      st.scaled = _scaled_count;
      if (cjIterations == 0 && !retry)
        _first_stats = st;
      if (_record_stats)
        _iteration_stats.push_back(st);
    }
    // the start-point policy may cache this iterate for the next solve
    if (middle)
      _start_policy.observeIteration(*this, st, cjIterations);

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
    // std::cout << ", " << cjIterations << ": " << _kappa << " "
    //       << FIXED(activeRobustChi2());
    ++cjIterations;
    postIteration(cjIterations);
    // std::cout << "*************** _kappa: " << _kappa << ", cjIterations: *"
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

    if (result == OptimizationAlgorithm::Fail) {
      return {cjIterations, true}; // optimize() returns 0, as before
    }

    bool converged = verifyConvergence(-1) &&
                     verifyEqFeasibility(_activeEdgesEq, -1.0) &&
                     verifyIneqFeasibility(_activeEdgesIneq, -1.0) &&
                     _kappa >= _kappa_final && auxiliariesValid();

    _last_converged = _last_converged || converged;
    stop = stop || converged;

    if (!stop) {
      computeKappa(kappa_flag);
    }
  }
  // std::cout <<
  // "***********************************************************************************************************************************************************************"
  //          << cjIterations << std::endl;
  if (verbose())
    std::cout << std::endl;
  return {cjIterations, false};
}

} // namespace cg2o
