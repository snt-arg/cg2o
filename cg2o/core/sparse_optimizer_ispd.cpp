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
  resetWarmStart();
  clearWarmStartShiftSequences();
  auto *oe = static_cast<g2o::OptimizableGraph::Edge *>(e);
  if (_edgeEqSet.count(oe) && !e->vertices().empty())
    _eq_multiplier_set.erase(
        static_cast<g2o::OptimizableGraph::Vertex *>(e->vertices().back()));
  edgeProcessingFunctionMap.erase(oe);
  getDualityGapFunctionMap.erase(oe);
  return SparseOptimizerForConstraints<SparseOptimizerISPD>::removeEdge(e);
}

bool SparseOptimizerISPD::removeVertex(g2o::HyperGraph::Vertex *v,
                                       bool detach) {
  if (!v || vertex(v->id()) != v)
    return false;
  resetWarmStart();
  clearWarmStartShiftSequences();
  auto *ov = static_cast<g2o::OptimizableGraph::Vertex *>(v);
  _snapshot_factories.erase(ov);
  _eq_multiplier_set.erase(ov);
  return SparseOptimizerForConstraints<SparseOptimizerISPD>::removeVertex(
      v, detach);
}

void SparseOptimizerISPD::clear() {
  resetWarmStart();
  clearWarmStartShiftSequences();
  _snapshot_factories.clear();
  _eq_multiplier_set.clear();
  edgeProcessingFunctionMap.clear();
  getDualityGapFunctionMap.clear();
  _last_converged = false;
  SparseOptimizerForConstraints<SparseOptimizerISPD>::clear();
}

void SparseOptimizerISPD::validateWarmStartSettings() const {
  auto positive = [](double x) { return warm_start::positiveFinite(x); };
  auto nonnegative = [](double x) { return x >= 0.0 && std::isfinite(x); };
  if (!positive(_kappa_initial) || !positive(_kappa_final) ||
      !positive(1.0 / _kappa_initial) || !positive(1.0 / _kappa_final) ||
      (_limit_kappa_final && _kappa_final < _kappa_initial) ||
      !positive(_kappa_update_factor) ||
      !positive(_slack_variable_initial_ineq) ||
      !positive(_lagrange_multiplier_initial_ineq))
    throw std::invalid_argument(
        "ISPD needs finite positive kappa and auxiliary initialization; "
        "kappa_final >= kappa_initial");
  if (_warm_start_strategy == 0)
    return;
  if (_warm_start_strategy < 0 || _warm_start_strategy > 3 ||
      _warm_start_shift_steps < 0 || _warm_start_middle_strategy < 0 ||
      _warm_start_middle_strategy > 3 || _warm_start_kappa_strategy < 0 ||
      _warm_start_kappa_strategy > 3 || _warm_start_reset_test < 0 ||
      _warm_start_reset_test > 2 || _warm_start_reset_mode < 0 ||
      _warm_start_reset_mode > 2 || !positive(_warm_start_middle_kappa) ||
      _warm_start_middle_iteration < 0 || !positive(_warm_start_middle_gap) ||
      !positive(_warm_start_kappa_factor) ||
      !positive(_warm_start_kappa_value) || !positive(_warm_start_kappa_gap) ||
      !nonnegative(_warm_start_push) || !nonnegative(_warm_start_slack_push) ||
      !nonnegative(_warm_start_multiplier_push) ||
      !nonnegative(_warm_start_reset_tolerance) ||
      !nonnegative(_warm_start_divergence_factor) ||
      !(_warm_start_center_lo >= 0.0 && _warm_start_center_lo <= 1.0 &&
        _warm_start_center_hi >= 1.0) ||
      !(_warm_start_middle_lo >= 0.0 && _warm_start_middle_lo <= 1.0 &&
        _warm_start_middle_hi >= 1.0) ||
      _warm_start_center_mode < 0 || _warm_start_center_mode > 3 ||
      _warm_start_slack_projection < 0 || _warm_start_slack_projection > 1)
    throw std::invalid_argument("invalid ISPD warm-start settings");
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

// Warm-start setters
void SparseOptimizerISPD::setWarmStartStrategy(int strategy) {
  if (strategy < 0 || strategy > 3) {
    throw std::invalid_argument("[Error] warm start strategy: 0 cold, 1 final, "
                                "2 middle, 3 current values");
  }
  _warm_start_strategy = strategy;
}

void SparseOptimizerISPD::setWarmStartShift(bool shift) {
  _warm_start_shift = shift;
}

void SparseOptimizerISPD::setWarmStartShiftSteps(int n) {
  if (n < 0) {
    throw std::invalid_argument("[Error] warm start shift steps must be >= 0");
  }
  _warm_start_shift_steps = n;
}

int SparseOptimizerISPD::warmStartShiftSteps() const {
  return _warm_start_shift_steps;
}

void SparseOptimizerISPD::setWarmStartPrimal(bool primal) {
  _warm_start_primal = primal;
}

void SparseOptimizerISPD::setWarmStartMiddleStrategy(int strategy) {
  if (strategy < 0 || strategy > 3) {
    throw std::invalid_argument("[Error] warm start middle strategy: 0 kappa, "
                                "1 duality-gap kappa, 2 iteration, "
                                "3 relative gap");
  }
  _warm_start_middle_strategy = strategy;
}

void SparseOptimizerISPD::setWarmStartMiddleKappa(double kappa) {
  _warm_start_middle_kappa = kappa;
}

void SparseOptimizerISPD::setWarmStartMiddleIteration(int iteration) {
  _warm_start_middle_iteration = iteration;
}

void SparseOptimizerISPD::setWarmStartKappaStrategy(int strategy) {
  if (strategy < 0 || strategy > 3) {
    throw std::invalid_argument("[Error] warm start kappa strategy: 0 stored, "
                                "1 duality gap, 2 constant, 3 relative gap");
  }
  _warm_start_kappa_strategy = strategy;
}

void SparseOptimizerISPD::setWarmStartKappaFactor(double factor) {
  _warm_start_kappa_factor = factor;
}

void SparseOptimizerISPD::setWarmStartKappaValue(double kappa) {
  _warm_start_kappa_value = kappa;
}

void SparseOptimizerISPD::setWarmStartPush(double push) {
  if (!(push >= 0.0) || !std::isfinite(push))
    throw std::invalid_argument("warm-start push must be finite and >= 0");
  _warm_start_push = push;
}

void SparseOptimizerISPD::setWarmStartSlackPush(double push) {
  if (!(push >= 0.0) || !std::isfinite(push))
    throw std::invalid_argument("slack push must be finite and >= 0");
  _warm_start_slack_push = push;
}

void SparseOptimizerISPD::setWarmStartMultiplierPush(double push) {
  if (!(push >= 0.0) || !std::isfinite(push))
    throw std::invalid_argument("multiplier push must be finite and >= 0");
  _warm_start_multiplier_push = push;
}

void SparseOptimizerISPD::setWarmStartSlackProjection(int strategy) {
  if (strategy < 0 || strategy > 1)
    throw std::invalid_argument("slack projection: 0 retain, 1 project to -g");
  _warm_start_slack_projection = strategy;
}

void SparseOptimizerISPD::setWarmStartCentering(double beta) {
  if (!std::isfinite(beta) || (beta != 0.0 && beta < 1.0)) {
    throw std::invalid_argument(
        "[Error] warm start centering: beta must be 0 (off) or >= 1");
  }
  if (beta == 0.0) {
    setWarmStartCenteringBand(0.0, std::numeric_limits<double>::infinity());
  } else {
    setWarmStartCenteringBand(1.0 / beta, beta);
  }
}

void SparseOptimizerISPD::setWarmStartCenteringBand(double lo, double hi) {
  if (!(lo >= 0.0) || !(hi >= 1.0) || lo > 1.0) {
    throw std::invalid_argument(
        "[Error] warm start centering band: need 0 <= lo <= 1 <= hi");
  }
  _warm_start_center_lo = lo;
  _warm_start_center_hi = hi;
}

void SparseOptimizerISPD::setWarmStartCenteringMode(int mode) {
  if (mode < 0 || mode > 3) {
    throw std::invalid_argument(
        "warm start centering mode: 0 ratio, "
        "1 legacy, 2 preserve slack, 3 preserve multiplier");
  }
  _warm_start_center_mode = mode;
}

bool SparseOptimizerISPD::centringActive() const {
  return _warm_start_center_lo > 0.0 ||
         _warm_start_center_hi < std::numeric_limits<double>::infinity();
}

void SparseOptimizerISPD::setWarmStartMiddleGap(double relativeGap) {
  _warm_start_middle_gap = relativeGap;
}

void SparseOptimizerISPD::setWarmStartMiddleNeighbourhood(double lo,
                                                          double hi) {
  if (!(lo >= 0.0) || !(hi >= 1.0) || lo > 1.0) {
    throw std::invalid_argument(
        "[Error] warm start middle neighbourhood: need 0 <= lo <= 1 <= hi");
  }
  _warm_start_middle_lo = lo;
  _warm_start_middle_hi = hi;
}

void SparseOptimizerISPD::setWarmStartKappaGap(double relativeGap) {
  _warm_start_kappa_gap = relativeGap;
}

bool SparseOptimizerISPD::inMiddleNeighbourhood(const PointStats &st) const {
  return st.rows == 0 || (st.minRatio >= _warm_start_middle_lo &&
                          st.maxRatio <= _warm_start_middle_hi);
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
  _acc_sum = 0.0;
  _acc_min = std::numeric_limits<double>::infinity();
  _acc_max = -std::numeric_limits<double>::infinity();
  _acc_rows = 0;
  for (auto &edge : _activeEdgesIneq) {
    executeEdgeProcessing(edge, 13);
  }
  st.rows = _acc_rows;
  st.kappa = _kappa;
  if (_acc_rows > 0) {
    st.mu = _acc_sum / _acc_rows;
    if (st.mu > 0.0) {
      st.minRatio = _acc_min / st.mu;
      st.maxRatio = _acc_max / st.mu;
      st.centrality =
          std::max(std::abs(st.minRatio - 1.0), std::abs(st.maxRatio - 1.0));
    }
  }
  if (withCost) {
    // Products and objective must describe the SAME iterate. Constraint
    // chi2() is zero; the robust objective also supports robust cost edges.
    computeActiveErrors();
    st.cost = 0.5 * activeRobustChi2();
    st.relativeGap = _acc_sum / (1.0 + std::abs(st.cost));
  }
  return st;
}

void SparseOptimizerISPD::setWarmStartFallback(bool fallback) {
  _warm_start_fallback = fallback;
}

void SparseOptimizerISPD::setWarmStartRetry(bool retry) {
  _warm_start_retry = retry;
}

void SparseOptimizerISPD::setWarmStartResetTolerance(double tolerance) {
  _warm_start_reset_tolerance = tolerance;
}

void SparseOptimizerISPD::setWarmStartDivergenceFactor(double factor) {
  _warm_start_divergence_factor = factor;
}

double SparseOptimizerISPD::warmStartDivergenceFactor() const {
  return _warm_start_divergence_factor;
}

bool SparseOptimizerISPD::lastRestored() const { return _last_restored; }

void SparseOptimizerISPD::setWarmStartResetTest(int test) {
  if (test < 0 || test > 2) {
    throw std::invalid_argument(
        "warm start reset test: 0 changed g, 1 residual with unit offset, 2 "
        "relative slack residual");
  }
  _warm_start_reset_test = test;
}

int SparseOptimizerISPD::warmStartResetTest() const {
  return _warm_start_reset_test;
}

void SparseOptimizerISPD::setWarmStartResetMode(int mode) {
  if (mode < 0 || mode > 2) {
    throw std::invalid_argument("[Error] warm start reset mode: 0 centred "
                                "rows, 1 cold rows, 2 cold duals");
  }
  _warm_start_reset_mode = mode;
}

int SparseOptimizerISPD::warmStartResetMode() const {
  return _warm_start_reset_mode;
}

// Warm-start getters
int SparseOptimizerISPD::warmStartStrategy() const {
  return _warm_start_strategy;
}
bool SparseOptimizerISPD::warmStartShift() const { return _warm_start_shift; }
bool SparseOptimizerISPD::warmStartPrimal() const { return _warm_start_primal; }
int SparseOptimizerISPD::warmStartMiddleStrategy() const {
  return _warm_start_middle_strategy;
}
double SparseOptimizerISPD::warmStartMiddleKappa() const {
  return _warm_start_middle_kappa;
}
int SparseOptimizerISPD::warmStartMiddleIteration() const {
  return _warm_start_middle_iteration;
}
int SparseOptimizerISPD::warmStartKappaStrategy() const {
  return _warm_start_kappa_strategy;
}
double SparseOptimizerISPD::warmStartKappaFactor() const {
  return _warm_start_kappa_factor;
}
double SparseOptimizerISPD::warmStartKappaValue() const {
  return _warm_start_kappa_value;
}
double SparseOptimizerISPD::warmStartPush() const { return _warm_start_push; }
double SparseOptimizerISPD::warmStartCentering() const {
  return centringActive() ? _warm_start_center_hi : 0.0;
}
bool SparseOptimizerISPD::warmStartFallback() const {
  return _warm_start_fallback;
}
bool SparseOptimizerISPD::warmStartRetry() const { return _warm_start_retry; }
double SparseOptimizerISPD::warmStartResetTolerance() const {
  return _warm_start_reset_tolerance;
}
bool SparseOptimizerISPD::lastRetried() const { return _last_retry; }
int SparseOptimizerISPD::lastResetRows() const { return _last_reset_rows; }

void SparseOptimizerISPD::clearWarmStartShiftSequences() {
  _shift_functions.clear();
  _shift_members.clear();
  _post_shift_callback = {};
}

void SparseOptimizerISPD::resetWarmStart() {
  discardVertexSnapshot();
  _warm_point_valid = false;
}

bool SparseOptimizerISPD::warmStartAvailable() const {
  return _warm_point_valid && (_last_converged || !_warm_start_fallback);
}

int SparseOptimizerISPD::lastWarmStartUsed() const { return _last_warm_used; }

double SparseOptimizerISPD::lastKappa() const { return _kappa; }

double SparseOptimizerISPD::lastStartKappa() const { return _last_start_kappa; }

bool SparseOptimizerISPD::lastConverged() const { return _last_converged; }

// Warm-start helpers

bool SparseOptimizerISPD::isEqMultiplierVertex(
    g2o::OptimizableGraph::Vertex *v) {
  return _eq_multiplier_set.count(v) > 0;
}

// Owned snapshots never touch g2o's push/pop stacks. Application backups and
// an algorithm's line-search backups may therefore coexist with warm data.
void SparseOptimizerISPD::discardSnapshot(
    std::vector<VertexSnapshotEntry> &list) {
  list.clear();
}

void SparseOptimizerISPD::restoreSnapshot(
    std::vector<VertexSnapshotEntry> &list, bool primal) {
  for (auto &entry : list) {
    auto *v = vertex(entry.id);
    if (v != entry.vertex || v->hessianIndex() < 0) {
      continue;
    }
    if (!v->fixed() && (entry.multiplier || primal)) {
      if (!entry.restore())
        throw std::runtime_error(
            "failed to restore warm-start vertex estimate");
    }
  }
  list.clear();
}

void SparseOptimizerISPD::pushSnapshot(std::vector<VertexSnapshotEntry> &list,
                                       bool primal) {
  discardSnapshot(list); // at most one entry per vertex and list
  _eq_multiplier_set.clear();
  _eq_multiplier_set.insert(_vEqLagrangeMultipliers.begin(),
                            _vEqLagrangeMultipliers.end());
  list.reserve(_activeVertices.size());
  for (auto *v : _activeVertices) {
    if (v->fixed()) {
      continue;
    }
    const bool multiplier = isEqMultiplierVertex(v);
    if (!multiplier && !primal) {
      continue; // primal not restored, so not stored
    }
    auto factory = _snapshot_factories.find(v);
    if (factory != _snapshot_factories.end()) {
      list.push_back({v, v->id(), multiplier, factory->second()});
    } else {
      std::vector<double> estimate;
      if (v->estimateDimension() <= 0 || !v->getEstimateData(estimate))
        throw std::invalid_argument("warm start needs typed addVertex() or "
                                    "registerWarmStartVertex<T>() for vertex " +
                                    std::to_string(v->id()));
      list.push_back(
          {v, v->id(), multiplier, [v, estimate = std::move(estimate)]() {
             return v->setEstimateData(estimate.data());
           }});
    }
  }
}

void SparseOptimizerISPD::discardVertexSnapshot() {
  discardSnapshot(_warm_snapshot);
  discardSnapshot(_start_snapshot);
}

void SparseOptimizerISPD::storeStartingPoint() {
  pushSnapshot(_start_snapshot, true);
  for (auto *edge : _activeEdgesIneq)
    executeEdgeProcessing(edge, 10);
  _start_point_kappa = _kappa;
}

void SparseOptimizerISPD::restoreStartingPoint() {
  restoreSnapshot(_start_snapshot, true);
  for (auto *edge : _activeEdgesIneq)
    executeEdgeProcessing(edge, 11);
  _kappa = _start_point_kappa;
  computeActiveErrors();
}

bool SparseOptimizerISPD::auxiliariesValid() {
  _aux_valid = true;
  for (auto *edge : _activeEdgesIneq)
    executeEdgeProcessing(edge, 15);
  return _aux_valid;
}

void SparseOptimizerISPD::storeWarmPoint(double kappa) {
  pushSnapshot(_warm_snapshot, _warm_start_primal);
  for (auto &edge : _activeEdgesIneq) {
    executeEdgeProcessing(edge, 4); // store lambda, s
  }
  _warm_point_kappa = kappa;
  _warm_point_stamp = _solve_counter;
  _warm_point_source = _warm_start_strategy;
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
  // _warm_start_lagrange_multiplier_eq_flag (shared with AL / BIPM) keeps the
  // equality multipliers of the previous solve
  if (forceEqualityReset || !_warm_start_lagrange_multiplier_eq_flag) {
    resetLagrangeMultiplierEq();
  }

  if (_activeEdgesIneq.empty()) {
    _kappa = _kappa_final;
  } else {
    _kappa = _kappa_initial;
  }

  for (auto &edge : _activeEdgesIneq) {
    // initialize the lagrange multipliers and slack variables
    executeEdgeProcessing(edge, 0);
  }

  if (_init_kappa_strategy == 1) {
    computeKappa(false);
  }
}

bool SparseOptimizerISPD::prepareStartingPoint() {
  _last_reset_rows = 0;
  const bool current = _warm_start_strategy == 3;
  const bool use = current || (_warm_start_strategy != 0 && _warm_point_valid &&
                               (_last_converged || !_warm_start_fallback));
  if (!use) {
    discardVertexSnapshot();
    coldStart();
    _last_warm_used = 0;
    return false;
  }

  auto shiftAll = [this]() {
    if (_warm_start_shift && _warm_start_shift_steps > 0) {
      for (auto &shift : _shift_functions) {
        shift(_warm_start_primal, _warm_start_shift_steps);
      }
      if (_post_shift_callback)
        _post_shift_callback();
    }
  };

  // 1. stored x (if primal warm start) and nu, stored lambda and s (edges
  //    without stored data are marked cold) - or, for strategy 3, the
  //    values currently in the graph; 2. MPC shift
  if (current) {
    discardVertexSnapshot();
    for (auto &edge : _activeEdgesIneq) {
      executeEdgeProcessing(edge, 14);
    }
  } else {
    restoreSnapshot(_warm_snapshot, _warm_start_primal);
    for (auto &edge : _activeEdgesIneq) {
      executeEdgeProcessing(edge, 5);
    }
  }
  shiftAll();
  const int source = current ? 3 : _warm_point_source;
  const double sourceKappa = (current && _kappa > 0.0)
                                 ? _kappa
                                 : (current ? _kappa_final : _warm_point_kappa);

  // drop the previous solve's starting point that lies below on the stacks
  discardSnapshot(_start_snapshot);

  if (_activeEdgesIneq.empty()) {
    _kappa = _kappa_final;
    _last_warm_used = source;
    return true;
  }

  auto clampKappa = [this, sourceKappa](double kappa) {
    if (!std::isfinite(kappa) || !(kappa > 0.0)) {
      kappa = warm_start::positiveFinite(sourceKappa) ? sourceKappa
                                                      : _kappa_initial;
    }
    kappa = std::max(kappa, _kappa_initial);
    if (_limit_kappa_final) {
      kappa = std::min(kappa, _kappa_final);
    }
    return kappa;
  };

  // 3. edges that were not active in the stored solve: cold initialisation
  //    (with the stored kappa, used by init_lagrange_strategy 1)
  _kappa = clampKappa(sourceKappa);
  for (auto &edge : _activeEdgesIneq) {
    executeEdgeProcessing(edge, 6);
  }

  // 4. change detector: rows whose stored (lambda, s) no longer fit the
  //    constraint at the warm primal (e.g. an obstacle appeared on the path)
  //    are re-initialised, or the whole solve starts with cold duals
  if (_warm_start_reset_tolerance > 0.0) {
    for (auto &edge : _activeEdgesIneq) {
      executeEdgeProcessing(edge, 9);
    }
    if (_last_reset_rows > 0 && _warm_start_reset_mode == 2) {
      // the problem changed too much for the stored duals: cold start of
      // nu, lambda, s and kappa from the warm primal
      coldStart(true);
      _last_warm_used = 0;
      return false;
    }
  }

  // 5. push s, lambda away from zero
  for (auto &edge : _activeEdgesIneq)
    executeEdgeProcessing(edge, 7);

  // 6. barrier parameter
  switch (_warm_start_kappa_strategy) {
  case 0:
    _kappa = clampKappa(_warm_start_kappa_factor * sourceKappa);
    break;
  case 1:
    _kappa = clampKappa(_warm_start_kappa_factor * dualityGapKappa());
    break;
  case 2:
    _kappa = clampKappa(_warm_start_kappa_value);
    break;
  case 3: {
    // mu0 = gap (1 + |cost|) / m at the warm primal
    computeActiveErrors();
    const double cost = 0.5 * activeRobustChi2();
    int m = 0;
    for (auto &edge : _activeEdgesIneq) {
      m += edge->dimension();
    }
    _kappa = clampKappa(m / (_warm_start_kappa_gap * (1.0 + std::abs(cost))));
    break;
  }
  default:
    throw std::runtime_error("[Error] Unknown warm start kappa strategy");
  }

  // 7. re-centre the pairs into the band around 1/kappa
  if (centringActive()) {
    for (auto &edge : _activeEdgesIneq) {
      executeEdgeProcessing(edge, 8);
    }
  }

  if (!auxiliariesValid()) {
    coldStart(true);
    _last_warm_used = 0;
    return false;
  }

  _last_warm_used = source;
  return true;
}

// solver functions

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
  // >= 1): before, they were then never computed and controller 3 applied
  // stale values.
  _step_size = ineq_backtracking_step; // for the prediction step
  _step_size_index = 0;

  for (auto &edge : _activeEdgesIneq) {
    // controller =1: compute the error prediction
    executeEdgeProcessing(edge, 1);
  }

  for (auto &edge : _activeEdgesIneq) {
    // controller =2: updates, and backtracking the lagrange multiplier and
    // slack variable (no backtracking when _aux_backtracking_step_min >= 1)
    executeEdgeProcessing(edge, 2);
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
    _last_retry = _last_restored = false;
    _last_warm_used = 0;
    resetWarmStart();
    return 0;
  }
  validateWarmStartSettings();
  if (_ivMap.size() == 0) {
    G2O_WARN("0 vertices to optimize, maybe forgot to call "
             "initializeOptimization()");
    return -1;
  }

  // starting point: cold (nu, lambda, s, kappa re-initialised) or warm
  // (stored point of the previous solve, optionally shifted), see
  // setWarmStartStrategy()
  ++_solve_counter;
  _last_retry = false;
  _last_restored = false;
  const bool warm = prepareStartingPoint(); // uses _last_converged
  _last_converged = false;
  _last_start_kappa = _kappa;
  _iteration_stats.clear();
  _start_stats = PointStats();
  _first_stats = PointStats();
  const bool wantStats = _warm_start_strategy != 0 || _record_stats;
  if (wantStats) {
    // cost needs the errors at the starting point (one extra evaluation)
    _start_stats = pointStats(_record_stats);
    _start_stats.iteration = -1;
  }
  const bool keepStart = _warm_start_strategy != 0;
  double startViolation = 0.0;
  if (keepStart) {
    // keep the starting point (x, nu) so that a failed solve can be undone
    storeStartingPoint();
    if (_warm_start_divergence_factor > 0.0) {
      startViolation = constraintViolation();
      _start_stats.violation = startViolation;
    }
  }

  // a solve broke down: linear-solver failure, or (not converged and) its
  // final constraint violation is far above the one it started from
  auto brokeDown = [&](bool failedSolve) {
    if (failedSolve || !auxiliariesValid()) {
      return true;
    }
    computeActiveErrors();
    if (!std::isfinite(activeRobustChi2()))
      return true;
    if (!keepStart || _last_converged || _warm_start_divergence_factor <= 0.0) {
      return false;
    }
    const double floor =
        std::max(startViolation, terminationCriteria().epsilon_constraint);
    return !(constraintViolation() <= _warm_start_divergence_factor * floor);
  };

  bool warm_stored = false; // warm point stored during this solve
  bool failed = false;
  int done = iterate(iterations, online, warm_stored, failed);
  _last_iteration_count = std::max(done, 0);
  if (done < 0) {
    if (keepStart) {
      restoreStartingPoint();
      _last_restored = true;
    }
    resetWarmStart();
    return -1;
  }

  bool broken = brokeDown(failed);
  if (broken && keepStart) {
    // put the starting point back instead of leaving the diverged iterate in
    // the vertices (a middle point stored by this solve is dropped)
    discardSnapshot(_warm_snapshot);
    restoreStartingPoint();
    _last_converged = false;
    warm_stored = false;
    _last_restored = true;
    if (warm && _warm_start_retry) {
      // and solve once more from it with a cold start of nu, lambda, s, kappa
      coldStart(true);
      storeStartingPoint();
      _last_warm_used = 0;
      _last_retry = true;
      failed = false;
      const int again = iterate(iterations, online, warm_stored, failed);
      _last_iteration_count += std::max(again, 0);
      if (again < 0) {
        restoreStartingPoint();
        _last_restored = true;
        resetWarmStart();
        return -1;
      }
      done += again;
      broken = brokeDown(failed);
      if (broken) {
        discardSnapshot(_warm_snapshot);
        restoreStartingPoint();
        _last_converged = false;
        warm_stored = false;
      } else {
        _last_restored = false;
      }
    }
  }
  if (failed || broken) {
    // the iterate is not usable: the next solve starts cold
    resetWarmStart();
    return failed ? 0 : done;
  }

  // warm start: store the final iterate (strategy 1), or the final one if
  // no middle iterate was stored (strategy 2)
  if (_warm_start_strategy == 1 ||
      (_warm_start_strategy == 2 && !warm_stored)) {
    storeWarmPoint(_kappa);
    warm_stored = true;
  }
  _warm_point_valid = _warm_start_strategy != 0 && warm_stored;
  discardSnapshot(_start_snapshot);
  if (!_warm_point_valid)
    discardSnapshot(_warm_snapshot);
  return done;
}

int SparseOptimizerISPD::iterate(int iterations, bool online, bool &warm_stored,
                                 bool &failed) {
  using namespace g2o;
  int cjIterations = 0;
  double cumTime = 0;
  bool ok = true;

  bool stop = false;
  failed = false;
  bool middle_centred = false; // a centred middle point was stored

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

    // statistics of the new iterate (first iteration: step sizes)
    const bool middle = ok && _warm_start_strategy == 2;
    const bool needStats = _record_stats || middle ||
                           (cjIterations == 0 && _warm_start_strategy != 0);
    PointStats st;
    if (needStats) {
      st = pointStats(_record_stats ||
                      (middle && _warm_start_middle_strategy == 3));
      st.iteration = cjIterations;
      st.attempt = _last_retry ? 1 : 0;
      st.ineqStep = _last_ineq_step;
      st.auxStep = _last_aux_step;
      st.appliedStep = _last_applied_step;
      st.scaled = _scaled_count;
      if (cjIterations == 0 && !_last_retry)
        _first_stats = st;
      if (_record_stats)
        _iteration_stats.push_back(st);
    }

    // warm start, strategy 2: keep the last iterate that qualifies as the
    // middle point (the first iterate if none does)
    if (middle && std::isfinite(st.mu) && st.mu > 0.0 && auxiliariesValid()) {
      bool qualifies = false;
      switch (_warm_start_middle_strategy) {
      case 0:
        qualifies = _kappa <= _warm_start_middle_kappa;
        break;
      case 1:
        qualifies =
            (st.mu > 0.0 ? 1.0 / st.mu : 0.0) <= _warm_start_middle_kappa;
        break;
      case 2:
        qualifies = cjIterations < _warm_start_middle_iteration;
        break;
      case 3:
        qualifies = st.relativeGap >= _warm_start_middle_gap;
        break;
      default:
        throw std::runtime_error("[Error] Unknown warm start middle strategy");
      }
      // the neighbourhood is a preference: a qualifying iterate inside it
      // replaces any stored one, a qualifying iterate outside it is stored
      // only while no centred one has been stored in this solve
      const bool centred = qualifies && inMiddleNeighbourhood(st);
      if (!warm_stored || centred || (qualifies && !middle_centred)) {
        storeWarmPoint(_kappa);
        warm_stored = true;
        middle_centred = middle_centred || centred;
      }
    }

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
      failed = true; // optimize() returns 0, as before
      return cjIterations;
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
  return cjIterations;
}

} // namespace cg2o
