#include "ispd_start_point.h"
#include "cg2o/core/sparse_optimizer_ispd.h"
#include "ispd_start_math.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace cg2o {

// ---- ISPDStartContext ----

void ISPDStartContext::coldStart() {
  _start.discardCopies();
  _start.initializeCold(false);
}

bool ISPDStartContext::storedPointUsable() const {
  return _start.storedPointUsable();
}

double ISPDStartContext::restoreStoredPoint() {
  return _start.restoreStoredPoint();
}

double ISPDStartContext::takeGraphPoint() { return _start.takeGraphPoint(); }

void ISPDStartContext::storeCurrentPoint() { _start.storeCurrentPoint(); }

bool ISPDStartContext::storedInThisSolve() const {
  return _start._stored_this_solve;
}

double ISPDStartContext::kappa() const { return _start.currentKappa(); }

bool ISPDStartContext::pairsValid() const { return _start.pairsValid(); }

const ISPDStartOptions &ISPDStartContext::options() const {
  return _start._options;
}

// ---- ISPDStartPoint ----

ISPDStartPoint::ISPDStartPoint(SparseOptimizerISPD &solver)
    : _solver(solver), _method(std::make_unique<ispd::Cold>()) {}

ISPDStartPoint::~ISPDStartPoint() = default;

void ISPDStartPoint::select(ISPDStartMethod method) {
  switch (method) {
  case ISPDStartMethod::Cold:
    use<ispd::Cold>();
    return;
  case ISPDStartMethod::PreviousFinal:
    use<ispd::PreviousFinal>();
    return;
  case ISPDStartMethod::PreviousMiddle:
    use<ispd::PreviousMiddle>();
    return;
  case ISPDStartMethod::CurrentGraph:
    use<ispd::CurrentGraph>();
    return;
  default:
    throw std::invalid_argument("ISPD start method: 0 cold, 1 previous "
                                "final, 2 previous middle, 3 current graph");
  }
}

void ISPDStartPoint::validate() const {
  const auto &o = _options;
  auto positive = [](double x) { return ispd::positiveFinite(x); };
  auto nonnegative = [](double x) { return x >= 0.0 && std::isfinite(x); };
  auto inRange = [](auto value, int last) {
    const int v = static_cast<int>(value);
    return v >= 0 && v <= last;
  };
  if (o.shift.steps < 0 || !inRange(o.kappa.rule, 3) ||
      !inRange(o.reset.test, 2) || !inRange(o.reset.mode, 2) ||
      !positive(o.kappa.factor) || !positive(o.kappa.value) ||
      !positive(o.kappa.gap) || !nonnegative(o.push.both) ||
      !nonnegative(o.push.slack) || !nonnegative(o.push.multiplier) ||
      !nonnegative(o.reset.tolerance) ||
      !nonnegative(o.recovery.divergenceFactor) ||
      !(o.centering.lo >= 0.0 && o.centering.lo <= 1.0 &&
        o.centering.hi >= 1.0) ||
      !inRange(o.centering.mode, 3))
    throw std::invalid_argument("invalid ISPD start-point options");
  _method.active().validate();
}

bool ISPDStartPoint::available() const {
  // the methods only read through a const context
  const ISPDStartContext context(const_cast<ISPDStartPoint &>(*this));
  return _method.active().available(context);
}

ISPDStartPoint::RowCache &
ISPDStartPoint::rowsOf(const g2o::OptimizableGraph::Edge *edge, int dimension) {
  auto it = _rows.find(edge);
  if (it == _rows.end())
    it = _rows.emplace(edge, RowCache(dimension)).first;
  return it->second;
}

bool ISPDStartPoint::isWarm(const g2o::OptimizableGraph::Edge *edge) const {
  auto it = _rows.find(edge);
  return it != _rows.end() && it->second.warm;
}

void ISPDStartPoint::shiftRows(const g2o::OptimizableGraph::Edge *to,
                               const g2o::OptimizableGraph::Edge *from,
                               int dimension) {
  auto source = _rows.find(from);
  if (source == _rows.end()) {
    // an edge without cached rows is cold
    auto target = _rows.find(to);
    if (target != _rows.end())
      target->second.warm = false;
    return;
  }
  const double *g = source->second.referenceConstraint();
  const bool warm = source->second.warm;
  RowCache &target = rowsOf(to, dimension);
  std::copy_n(g, dimension, target.referenceConstraint());
  target.warm = warm;
}

void ISPDStartPoint::clearShiftSequences() {
  _shift_functions.clear();
  _shift_members.clear();
  _post_shift_callback = {};
}

void ISPDStartPoint::discardCopies() {
  _vertices.discard(VertexEstimateCache::Cached);
  _vertices.discard(VertexEstimateCache::Rollback);
}

void ISPDStartPoint::invalidate() {
  discardCopies();
  _warm_point_valid = false;
}

void ISPDStartPoint::forgetVertex(const g2o::OptimizableGraph::Vertex *v) {
  invalidate();
  clearShiftSequences();
  _vertices.forgetVertex(v);
}

void ISPDStartPoint::forgetEdge(const g2o::OptimizableGraph::Edge *e) {
  invalidate();
  clearShiftSequences();
  _rows.erase(e);
}

void ISPDStartPoint::clear() {
  invalidate();
  clearShiftSequences();
  _vertices.clear();
  _rows.clear();
}

void ISPDStartPoint::skipSolve() {
  _report = ISPDStartReport();
  invalidate();
}

// ---- the solver's cold initialization rules ----

void ISPDStartPoint::initializeSlack(IneqBlock &block) {
  const double initial = _solver._slack_variable_initial_ineq;
  switch (_solver._init_slack_strategy) {
  case 0:   // max(-g(x0), initial)
  case 1: { // max(|g(x0)|, initial)
    block.edge->computeError();
    const double *g = block.edge->errorData();
    for (int i = 0; i < block.dimension; ++i) {
      const double candidate =
          _solver._init_slack_strategy == 0 ? -g[i] : std::abs(g[i]);
      block.slack[i] = candidate > initial ? candidate : initial;
    }
    break;
  }
  case 2: // constant
    std::fill_n(block.slack, block.dimension, initial);
    break;
  default:
    throw std::runtime_error("[Error] Unknown slack initialization strategy. "
                             "0: max(-g(x0),_slack_variable_initial_ineq)"
                             "1: max(abs(g(x0)),_slack_variable_initial_ineq)"
                             "2: constant  =  _slack_variable_initial_ineq");
  }
}

void ISPDStartPoint::initializeColdRows(IneqBlock &block) {
  initializeSlack(block);
  const double initial = _solver._lagrange_multiplier_initial_ineq;
  switch (_solver._init_lagrange_strategy) {
  case 0:
    std::fill_n(block.multiplier, block.dimension, initial);
    break;
  case 1:
    for (int i = 0; i < block.dimension; ++i)
      block.multiplier[i] =
          std::max(1.0 / (_solver._kappa * block.slack[i]), initial);
    break;
  }
}

void ISPDStartPoint::initializeNewRows(IneqBlock &block) {
  // the multiplier rule may need a kappa, which is set when a solve starts
  initializeSlack(block);
  std::fill_n(block.multiplier, block.dimension,
              _solver._lagrange_multiplier_initial_ineq);
}

void ISPDStartPoint::initializeCold(bool forceEqualityReset) {
  // keepEqualityMultipliers keeps the equality multipliers of the previous
  // solve
  if (forceEqualityReset || !_options.keepEqualityMultipliers)
    _solver.resetLagrangeMultiplierEq();
  _solver._kappa = _solver._activeEdgesIneq.empty() ? _solver._kappa_final
                                                    : _solver._kappa_initial;
  _solver.forEachActiveIneqBlock(
      [&](IneqBlock &block) { initializeColdRows(block); });
  if (_solver._init_kappa_strategy == 1)
    _solver.computeKappa(false);
}

// ---- stored and rollback points ----

bool ISPDStartPoint::storedPointUsable() const {
  return _warm_point_valid && _warm_point_method == _method.type() &&
         (_solver._last_converged || !_options.recovery.fallback);
}

void ISPDStartPoint::storeCurrentPoint() {
  _vertices.store(VertexEstimateCache::Cached, _solver._activeVertices,
                  _solver._vEqLagrangeMultipliers, _options.warmPrimal);
  _solver.forEachActiveIneqBlock([&](IneqBlock &block) {
    block.edge->computeError();
    RowCache &rows = rowsOf(block.edge, block.dimension);
    std::copy_n(block.slack, block.dimension, rows.cachedSlack());
    std::copy_n(block.multiplier, block.dimension, rows.cachedMultiplier());
    std::copy_n(block.edge->errorData(), block.dimension,
                rows.storedConstraint());
    rows.cachedStamp = _solve_counter;
  });
  _warm_point_kappa = _solver._kappa;
  _warm_point_stamp = _solve_counter;
  _warm_point_method = _method.type();
  _stored_this_solve = true;
}

double ISPDStartPoint::restoreStoredPoint() {
  // stored x (if primal warm start) and nu, stored lambda and s; edges
  // without stored rows are marked cold
  _vertices.restore(VertexEstimateCache::Cached, _solver, _options.warmPrimal);
  _solver.forEachActiveIneqBlock([&](IneqBlock &block) {
    auto it = _rows.find(block.edge);
    if (it == _rows.end())
      return; // never stored: cold
    RowCache &rows = it->second;
    rows.warm = rows.cachedStamp == _warm_point_stamp;
    if (!rows.warm)
      return;
    std::copy_n(rows.cachedSlack(), block.dimension, block.slack);
    std::copy_n(rows.cachedMultiplier(), block.dimension, block.multiplier);
    std::copy_n(rows.storedConstraint(), block.dimension,
                rows.referenceConstraint());
  });
  return _warm_point_kappa;
}

double ISPDStartPoint::takeGraphPoint() {
  discardCopies();
  _solver.forEachActiveIneqBlock([&](IneqBlock &block) {
    block.edge->computeError();
    RowCache &rows = rowsOf(block.edge, block.dimension);
    std::copy_n(block.edge->errorData(), block.dimension,
                rows.referenceConstraint());
    rows.warm = true;
  });
  return _solver._kappa > 0.0 ? _solver._kappa : _solver._kappa_final;
}

double ISPDStartPoint::currentKappa() const { return _solver._kappa; }

bool ISPDStartPoint::pairsValid() const { return _solver.auxiliariesValid(); }

void ISPDStartPoint::storeRollback() {
  _vertices.store(VertexEstimateCache::Rollback, _solver._activeVertices,
                  _solver._vEqLagrangeMultipliers, true);
  _solver.forEachActiveIneqBlock([&](IneqBlock &block) {
    RowCache &rows = rowsOf(block.edge, block.dimension);
    std::copy_n(block.slack, block.dimension, rows.rollbackSlack());
    std::copy_n(block.multiplier, block.dimension, rows.rollbackMultiplier());
    rows.rollbackStamp = _solve_counter;
  });
  _rollback_kappa = _solver._kappa;
}

void ISPDStartPoint::restoreRollback() {
  _vertices.restore(VertexEstimateCache::Rollback, _solver, true);
  _solver.forEachActiveIneqBlock([&](IneqBlock &block) {
    auto it = _rows.find(block.edge);
    if (it == _rows.end() || it->second.rollbackStamp != _solve_counter)
      return;
    std::copy_n(it->second.rollbackSlack(), block.dimension, block.slack);
    std::copy_n(it->second.rollbackMultiplier(), block.dimension,
                block.multiplier);
  });
  _solver._kappa = _rollback_kappa;
  _solver.computeActiveErrors();
}

// ---- one optimize() ----

bool ISPDStartPoint::prepareStartingPoint() {
  // 1. the method's point: stored x (with warmPrimal), nu, lambda and s, the
  //    values in the graph, or a cold start that is used as it is
  ISPDStartContext context(*this);
  const std::optional<double> loaded = _method.active().load(context);
  if (!loaded) {
    _report.used = ISPDStartMethod::Cold;
    return false;
  }
  // 2. MPC shift
  if (_options.shift.enabled && _options.shift.steps > 0) {
    for (auto &shift : _shift_functions)
      shift(_options.warmPrimal, _options.shift.steps);
    if (_post_shift_callback)
      _post_shift_callback();
  }
  const ISPDStartMethod source = _method.active().method();
  const double sourceKappa = *loaded;

  // drop the previous solve's rollback point
  _vertices.discard(VertexEstimateCache::Rollback);

  if (_solver._activeEdgesIneq.empty()) {
    _solver._kappa = _solver._kappa_final;
    _report.used = source;
    return true;
  }

  auto clampKappa = [this, sourceKappa](double kappa) {
    if (!std::isfinite(kappa) || !(kappa > 0.0)) {
      kappa = ispd::positiveFinite(sourceKappa) ? sourceKappa
                                                : _solver._kappa_initial;
    }
    kappa = std::max(kappa, _solver._kappa_initial);
    if (_solver._limit_kappa_final) {
      kappa = std::min(kappa, _solver._kappa_final);
    }
    return kappa;
  };

  // 3. edges without warm rows: cold initialisation (with the stored kappa,
  //    used by init_lagrange_strategy 1)
  _solver._kappa = clampKappa(sourceKappa);
  _solver.forEachActiveIneqBlock([&](IneqBlock &block) {
    if (!isWarm(block.edge))
      initializeColdRows(block);
  });

  // 4. change detector: rows whose stored (lambda, s) no longer fit the
  //    constraint at the warm primal (e.g. an obstacle appeared on the path)
  //    are re-initialised, or the whole solve starts with cold duals
  const auto &reset = _options.reset;
  if (reset.tolerance > 0.0) {
    _solver.forEachActiveIneqBlock([&](IneqBlock &block) {
      auto it = _rows.find(block.edge);
      if (it == _rows.end() || !it->second.warm)
        return; // already initialised cold
      block.edge->computeError();
      const double *g = block.edge->errorData();
      const double *stored = it->second.referenceConstraint();
      for (int i = 0; i < block.dimension; ++i)
        if (ispd::resetChangedRow(block.slack[i], block.multiplier[i], g[i],
                                  stored[i], _solver._kappa,
                                  static_cast<int>(reset.test), reset.tolerance,
                                  static_cast<int>(reset.mode),
                                  _solver._slack_variable_initial_ineq,
                                  _solver._lagrange_multiplier_initial_ineq))
          ++_report.resetRows;
    });
    if (_report.resetRows > 0 && reset.mode == ResetMode::ColdDuals) {
      // the problem changed too much for the stored duals: cold start of
      // nu, lambda, s and kappa from the warm primal
      initializeCold(true);
      _report.used = ISPDStartMethod::Cold;
      return false;
    }
  }

  // 5. repair invalid pairs, push s, lambda away from zero
  const double slackFloor = _options.push.slackFloor();
  const double multiplierFloor = _options.push.multiplierFloor();
  _solver.forEachActiveIneqBlock([&](IneqBlock &block) {
    block.edge->computeError();
    const double *g = block.edge->errorData();
    for (int i = 0; i < block.dimension; ++i)
      if (ispd::pushRow(block.slack[i], block.multiplier[i], g[i],
                        _solver._kappa, slackFloor, multiplierFloor,
                        _solver._slack_variable_initial_ineq,
                        _solver._lagrange_multiplier_initial_ineq,
                        _options.push.projectSlack))
        ++_report.resetRows;
  });

  // 6. barrier parameter
  const auto &k = _options.kappa;
  switch (k.rule) {
  case WarmKappaRule::PreviousTimesFactor:
    _solver._kappa = clampKappa(k.factor * sourceKappa);
    break;
  case WarmKappaRule::DualityGap:
    _solver._kappa = clampKappa(k.factor * _solver.dualityGapKappa());
    break;
  case WarmKappaRule::Constant:
    _solver._kappa = clampKappa(k.value);
    break;
  case WarmKappaRule::RelativeGap: {
    // mu0 = gap (1 + |cost|) / m at the warm primal
    _solver.computeActiveErrors();
    const double cost = 0.5 * _solver.activeRobustChi2();
    int m = 0;
    for (auto &edge : _solver._activeEdgesIneq) {
      m += edge->dimension();
    }
    _solver._kappa = clampKappa(m / (k.gap * (1.0 + std::abs(cost))));
    break;
  }
  default:
    throw std::runtime_error("[Error] Unknown warm start kappa strategy");
  }

  // 7. re-centre the pairs into the band around 1/kappa; a pair that cannot
  //    be centred fails the validity check below
  const auto &c = _options.centering;
  if (c.active()) {
    const double mu = 1.0 / _solver._kappa;
    _solver.forEachActiveIneqBlock([&](IneqBlock &block) {
      for (int i = 0; i < block.dimension; ++i)
        ispd::centerPair(block.slack[i], block.multiplier[i], mu, c.lo, c.hi,
                         static_cast<int>(c.mode), slackFloor, multiplierFloor);
    });
  }

  if (!_solver.auxiliariesValid()) {
    initializeCold(true);
    _report.used = ISPDStartMethod::Cold;
    return false;
  }

  _report.used = source;
  return true;
}

void ISPDStartPoint::prepare() {
  ++_solve_counter;
  _report = ISPDStartReport();
  _retrying = false;
  _stored_this_solve = false;
  _start_violation = 0.0;
  _method.active().beginAttempt();

  _started_warm = prepareStartingPoint(); // uses the last convergence
  _keep_rollback = _method.active().warm();
  _report.kappa = _solver._kappa;
  double startViolation = std::numeric_limits<double>::quiet_NaN();
  if (_keep_rollback) {
    // keep the starting point (x, nu, s, lambda, kappa) so that a failed
    // solve can be undone
    storeRollback();
    if (_options.recovery.divergenceFactor > 0.0) {
      _start_violation = _solver.constraintViolation();
      startViolation = _start_violation;
    }
  }
  if (_method.active().warm() || _solver._record_stats) {
    // the cost needs the errors at the starting point (one extra evaluation)
    _report.start = _solver.pointStats(_solver._record_stats);
    _report.start.iteration = -1;
  }
  _report.start.violation = startViolation;
}

ISPDStartPoint::StatsRequest ISPDStartPoint::statsRequest(int iteration,
                                                          bool accepted) const {
  // the first iteration of a warm method is reported; a method may observe
  // every accepted iterate
  const IterateObservation observation = _method.active().observation();
  const bool observe = accepted && observation != IterateObservation::None;
  return {observe || (iteration == 0 && _method.active().warm()),
          observe && observation == IterateObservation::ComplementarityAndCost};
}

void ISPDStartPoint::afterIteration(const ISPDPointStats &stats,
                                    bool accepted) {
  if (stats.iteration == 0 && stats.attempt == 0)
    _report.firstIteration = stats;
  if (accepted && _method.active().observation() != IterateObservation::None) {
    ISPDStartContext context(*this);
    _method.active().observe(context, stats);
  }
}

// A solve broke down: linear-solver failure, invalid auxiliaries, a
// non-finite objective, or (not converged and) a final constraint violation
// far above the one it started from.
bool ISPDStartPoint::brokeDown(bool solverFailed) {
  if (solverFailed || !_solver.auxiliariesValid()) {
    return true;
  }
  _solver.computeActiveErrors();
  if (!std::isfinite(_solver.activeRobustChi2()))
    return true;
  const double factor = _options.recovery.divergenceFactor;
  if (!_keep_rollback || _solver._last_converged || factor <= 0.0) {
    return false;
  }
  const double floor = std::max(
      _start_violation, _solver.terminationCriteria().epsilon_constraint);
  return !(_solver.constraintViolation() <= factor * floor);
}

ISPDStartPoint::FinishAction
ISPDStartPoint::finish(const SolveOutcome &outcome) {
  if (outcome.iterations < 0) {
    if (_keep_rollback) {
      restoreRollback();
      _report.restored = true;
    }
    invalidate();
    return FinishAction::Error;
  }

  const bool broken = brokeDown(outcome.solverFailed);
  if (broken && _keep_rollback) {
    // put the starting point back instead of leaving the diverged iterate in
    // the vertices (a point stored by this solve is dropped)
    _vertices.discard(VertexEstimateCache::Cached);
    restoreRollback();
    _solver._last_converged = false;
    _stored_this_solve = false;
    _report.restored = true;
    if (!_retrying && _started_warm && _options.recovery.retry)
      return FinishAction::RetryCold;
  } else if (_retrying) {
    _report.restored = false; // the cold retry recovered
  }
  if (outcome.solverFailed || broken) {
    // the iterate is not usable: the next solve starts cold
    invalidate();
    return FinishAction::Reject;
  }

  // the method keeps the point of the next solve
  ISPDStartContext context(*this);
  _method.active().retain(context);
  _warm_point_valid = _stored_this_solve;
  _vertices.discard(VertexEstimateCache::Rollback);
  if (!_warm_point_valid)
    _vertices.discard(VertexEstimateCache::Cached);
  return FinishAction::Accept;
}

void ISPDStartPoint::prepareColdRetry() {
  // solve once more from the restored point with a cold start of nu,
  // lambda, s and kappa
  initializeCold(true);
  storeRollback();
  _report.used = ISPDStartMethod::Cold;
  _report.retried = true;
  _retrying = true;
  _method.active().beginAttempt();
}

} // namespace cg2o
