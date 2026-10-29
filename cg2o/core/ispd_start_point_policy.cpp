#include "ispd_start_point_policy.h"
#include "sparse_optimizer_ispd.h"
#include "warm_start_math.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace cg2o {

void ISPDStartPointPolicy::validate() const {
  const auto &o = _options;
  auto positive = [](double x) { return warm_start::positiveFinite(x); };
  auto nonnegative = [](double x) { return x >= 0.0 && std::isfinite(x); };
  auto inRange = [](auto value, int last) {
    const int v = static_cast<int>(value);
    return v >= 0 && v <= last;
  };
  if (!inRange(o.source, 3) || o.shift.steps < 0 ||
      !inRange(o.middle.rule, 3) || !inRange(o.kappa.rule, 3) ||
      !inRange(o.reset.test, 2) || !inRange(o.reset.mode, 2) ||
      !positive(o.middle.kappa) || o.middle.iteration < 0 ||
      !positive(o.middle.gap) || !positive(o.kappa.factor) ||
      !positive(o.kappa.value) || !positive(o.kappa.gap) ||
      !nonnegative(o.push.both) || !nonnegative(o.push.slack) ||
      !nonnegative(o.push.multiplier) || !nonnegative(o.reset.tolerance) ||
      !nonnegative(o.recovery.divergenceFactor) ||
      !(o.centering.lo >= 0.0 && o.centering.lo <= 1.0 &&
        o.centering.hi >= 1.0) ||
      !(o.middle.neighbourhoodLo >= 0.0 && o.middle.neighbourhoodLo <= 1.0 &&
        o.middle.neighbourhoodHi >= 1.0) ||
      !inRange(o.centering.mode, 3))
    throw std::invalid_argument("invalid ISPD warm-start settings");
}

ISPDStartPointPolicy::RowCache &
ISPDStartPointPolicy::rowsOf(const g2o::OptimizableGraph::Edge *edge,
                             int dimension) {
  auto it = _rows.find(edge);
  if (it == _rows.end())
    it = _rows.emplace(edge, RowCache(dimension)).first;
  return it->second;
}

bool ISPDStartPointPolicy::isWarm(
    const g2o::OptimizableGraph::Edge *edge) const {
  auto it = _rows.find(edge);
  return it != _rows.end() && it->second.warm;
}

void ISPDStartPointPolicy::shiftRows(const g2o::OptimizableGraph::Edge *to,
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

void ISPDStartPointPolicy::clearShiftSequences() {
  _shift_functions.clear();
  _shift_members.clear();
  _post_shift_callback = {};
}

void ISPDStartPointPolicy::invalidate() {
  _vertices.discard(VertexEstimateCache::Cached);
  _vertices.discard(VertexEstimateCache::Rollback);
  _warm_point_valid = false;
}

void ISPDStartPointPolicy::forgetVertex(
    const g2o::OptimizableGraph::Vertex *v) {
  invalidate();
  clearShiftSequences();
  _vertices.forgetVertex(v);
}

void ISPDStartPointPolicy::forgetEdge(const g2o::OptimizableGraph::Edge *e) {
  invalidate();
  clearShiftSequences();
  _rows.erase(e);
}

void ISPDStartPointPolicy::clear() {
  invalidate();
  clearShiftSequences();
  _vertices.clear();
  _rows.clear();
}

void ISPDStartPointPolicy::skipSolve() {
  _last_retry = _last_restored = false;
  _last_source = 0;
  invalidate();
}

void ISPDStartPointPolicy::storeCachedPoint(SparseOptimizerISPD &solver,
                                            double kappa) {
  _vertices.store(VertexEstimateCache::Cached, solver._activeVertices,
                  solver._vEqLagrangeMultipliers, _options.warmPrimal);
  solver.forEachActiveIneqBlock([&](IneqBlock &block) {
    block.edge->computeError();
    RowCache &rows = rowsOf(block.edge, block.dimension);
    std::copy_n(block.slack, block.dimension, rows.cachedSlack());
    std::copy_n(block.multiplier, block.dimension, rows.cachedMultiplier());
    std::copy_n(block.edge->errorData(), block.dimension,
                rows.storedConstraint());
    rows.cachedStamp = _solve_counter;
  });
  _warm_point_kappa = kappa;
  _warm_point_stamp = _solve_counter;
  _warm_point_source = static_cast<int>(_options.source);
}

void ISPDStartPointPolicy::storeRollback(SparseOptimizerISPD &solver) {
  _vertices.store(VertexEstimateCache::Rollback, solver._activeVertices,
                  solver._vEqLagrangeMultipliers, true);
  solver.forEachActiveIneqBlock([&](IneqBlock &block) {
    RowCache &rows = rowsOf(block.edge, block.dimension);
    std::copy_n(block.slack, block.dimension, rows.rollbackSlack());
    std::copy_n(block.multiplier, block.dimension, rows.rollbackMultiplier());
    rows.rollbackStamp = _solve_counter;
  });
  _rollback_kappa = solver._kappa;
}

void ISPDStartPointPolicy::restoreRollback(SparseOptimizerISPD &solver) {
  _vertices.restore(VertexEstimateCache::Rollback, solver, true);
  solver.forEachActiveIneqBlock([&](IneqBlock &block) {
    auto it = _rows.find(block.edge);
    if (it == _rows.end() || it->second.rollbackStamp != _solve_counter)
      return;
    std::copy_n(it->second.rollbackSlack(), block.dimension, block.slack);
    std::copy_n(it->second.rollbackMultiplier(), block.dimension,
                block.multiplier);
  });
  solver._kappa = _rollback_kappa;
  solver.computeActiveErrors();
}

bool ISPDStartPointPolicy::prepareStartingPoint(SparseOptimizerISPD &solver) {
  _last_reset_rows = 0;
  const bool current = _options.source == ISPDStartSource::CurrentGraph;
  if (!current && !available(solver._last_converged)) {
    _vertices.discard(VertexEstimateCache::Cached);
    _vertices.discard(VertexEstimateCache::Rollback);
    solver.coldStart();
    _last_source = 0;
    return false;
  }

  // 1. stored x (if primal warm start) and nu, stored lambda and s (edges
  //    without stored rows are marked cold) - or, for CurrentGraph, the
  //    values currently in the graph; 2. MPC shift
  if (current) {
    _vertices.discard(VertexEstimateCache::Cached);
    _vertices.discard(VertexEstimateCache::Rollback);
    solver.forEachActiveIneqBlock([&](IneqBlock &block) {
      block.edge->computeError();
      RowCache &rows = rowsOf(block.edge, block.dimension);
      std::copy_n(block.edge->errorData(), block.dimension,
                  rows.referenceConstraint());
      rows.warm = true;
    });
  } else {
    _vertices.restore(VertexEstimateCache::Cached, solver, _options.warmPrimal);
    solver.forEachActiveIneqBlock([&](IneqBlock &block) {
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
  }
  if (_options.shift.enabled && _options.shift.steps > 0) {
    for (auto &shift : _shift_functions)
      shift(_options.warmPrimal, _options.shift.steps);
    if (_post_shift_callback)
      _post_shift_callback();
  }
  const int source = current ? 3 : _warm_point_source;
  const double sourceKappa =
      (current && solver._kappa > 0.0)
          ? solver._kappa
          : (current ? solver._kappa_final : _warm_point_kappa);

  // drop the previous solve's rollback point
  _vertices.discard(VertexEstimateCache::Rollback);

  if (solver._activeEdgesIneq.empty()) {
    solver._kappa = solver._kappa_final;
    _last_source = source;
    return true;
  }

  auto clampKappa = [&solver, sourceKappa](double kappa) {
    if (!std::isfinite(kappa) || !(kappa > 0.0)) {
      kappa = warm_start::positiveFinite(sourceKappa) ? sourceKappa
                                                      : solver._kappa_initial;
    }
    kappa = std::max(kappa, solver._kappa_initial);
    if (solver._limit_kappa_final) {
      kappa = std::min(kappa, solver._kappa_final);
    }
    return kappa;
  };

  // 3. edges without warm rows: cold initialisation (with the stored kappa,
  //    used by init_lagrange_strategy 1)
  solver._kappa = clampKappa(sourceKappa);
  solver.forEachActiveIneqBlock([&](IneqBlock &block) {
    if (!isWarm(block.edge))
      solver.executeEdgeProcessing(block.edge,
                                   IneqStateOperation::ColdInitialize);
  });

  // 4. change detector: rows whose stored (lambda, s) no longer fit the
  //    constraint at the warm primal (e.g. an obstacle appeared on the path)
  //    are re-initialised, or the whole solve starts with cold duals
  const auto &reset = _options.reset;
  if (reset.tolerance > 0.0) {
    solver.forEachActiveIneqBlock([&](IneqBlock &block) {
      auto it = _rows.find(block.edge);
      if (it == _rows.end() || !it->second.warm)
        return; // already initialised cold
      block.edge->computeError();
      const double *g = block.edge->errorData();
      const double *stored = it->second.referenceConstraint();
      for (int i = 0; i < block.dimension; ++i)
        if (warm_start::resetChangedRow(
                block.slack[i], block.multiplier[i], g[i], stored[i],
                solver._kappa, static_cast<int>(reset.test), reset.tolerance,
                static_cast<int>(reset.mode),
                solver._slack_variable_initial_ineq,
                solver._lagrange_multiplier_initial_ineq))
          ++_last_reset_rows;
    });
    if (_last_reset_rows > 0 && reset.mode == ResetMode::ColdDuals) {
      // the problem changed too much for the stored duals: cold start of
      // nu, lambda, s and kappa from the warm primal
      solver.coldStart(true);
      _last_source = 0;
      return false;
    }
  }

  // 5. repair invalid pairs, push s, lambda away from zero
  const double slackFloor = _options.push.slackFloor();
  const double multiplierFloor = _options.push.multiplierFloor();
  solver.forEachActiveIneqBlock([&](IneqBlock &block) {
    block.edge->computeError();
    const double *g = block.edge->errorData();
    for (int i = 0; i < block.dimension; ++i)
      if (warm_start::pushRow(block.slack[i], block.multiplier[i], g[i],
                              solver._kappa, slackFloor, multiplierFloor,
                              solver._slack_variable_initial_ineq,
                              solver._lagrange_multiplier_initial_ineq,
                              _options.push.projectSlack))
        ++_last_reset_rows;
  });

  // 6. barrier parameter
  const auto &k = _options.kappa;
  switch (k.rule) {
  case WarmKappaRule::PreviousTimesFactor:
    solver._kappa = clampKappa(k.factor * sourceKappa);
    break;
  case WarmKappaRule::DualityGap:
    solver._kappa = clampKappa(k.factor * solver.dualityGapKappa());
    break;
  case WarmKappaRule::Constant:
    solver._kappa = clampKappa(k.value);
    break;
  case WarmKappaRule::RelativeGap: {
    // mu0 = gap (1 + |cost|) / m at the warm primal
    solver.computeActiveErrors();
    const double cost = 0.5 * solver.activeRobustChi2();
    int m = 0;
    for (auto &edge : solver._activeEdgesIneq) {
      m += edge->dimension();
    }
    solver._kappa = clampKappa(m / (k.gap * (1.0 + std::abs(cost))));
    break;
  }
  default:
    throw std::runtime_error("[Error] Unknown warm start kappa strategy");
  }

  // 7. re-centre the pairs into the band around 1/kappa; a pair that cannot
  //    be centred fails the validity check below
  const auto &c = _options.centering;
  if (c.active()) {
    const double mu = 1.0 / solver._kappa;
    solver.forEachActiveIneqBlock([&](IneqBlock &block) {
      for (int i = 0; i < block.dimension; ++i)
        warm_start::centerPair(block.slack[i], block.multiplier[i], mu, c.lo,
                               c.hi, static_cast<int>(c.mode), slackFloor,
                               multiplierFloor);
    });
  }

  if (!solver.auxiliariesValid()) {
    solver.coldStart(true);
    _last_source = 0;
    return false;
  }

  _last_source = source;
  return true;
}

StartResult ISPDStartPointPolicy::prepare(SparseOptimizerISPD &solver) {
  ++_solve_counter;
  _last_retry = false;
  _last_restored = false;
  _retrying = false;
  _stored_this_solve = false;
  _middle_centred = false;
  _start_violation = 0.0;

  StartResult result;
  result.warm = prepareStartingPoint(solver); // uses the last convergence
  _started_warm = result.warm;
  _keep_rollback = _options.source != ISPDStartSource::Cold;
  result.keepRollback = _keep_rollback;
  if (_keep_rollback) {
    // keep the starting point (x, nu, s, lambda, kappa) so that a failed
    // solve can be undone
    storeRollback(solver);
    if (_options.recovery.divergenceFactor > 0.0) {
      _start_violation = solver.constraintViolation();
      result.startViolation = _start_violation;
    }
  }
  return result;
}

bool ISPDStartPointPolicy::inMiddleNeighbourhood(
    const ISPDPointStats &stats) const {
  return stats.rows == 0 ||
         (stats.minRatio >= _options.middle.neighbourhoodLo &&
          stats.maxRatio <= _options.middle.neighbourhoodHi);
}

void ISPDStartPointPolicy::observeIteration(SparseOptimizerISPD &solver,
                                            const ISPDPointStats &stats,
                                            int iteration) {
  // keep the last iterate that qualifies as the middle point (the first
  // iterate if none does)
  if (!(std::isfinite(stats.mu) && stats.mu > 0.0 && solver.auxiliariesValid()))
    return;
  const auto &m = _options.middle;
  bool qualifies = false;
  switch (m.rule) {
  case MiddlePointRule::BarrierKappa:
    qualifies = solver._kappa <= m.kappa;
    break;
  case MiddlePointRule::DualityGapKappa:
    qualifies = (stats.mu > 0.0 ? 1.0 / stats.mu : 0.0) <= m.kappa;
    break;
  case MiddlePointRule::Iteration:
    qualifies = iteration < m.iteration;
    break;
  case MiddlePointRule::RelativeGap:
    qualifies = stats.relativeGap >= m.gap;
    break;
  default:
    throw std::runtime_error("[Error] Unknown warm start middle strategy");
  }
  // the neighbourhood is a preference: a qualifying iterate inside it
  // replaces any stored one, a qualifying iterate outside it is stored
  // only while no centred one has been stored in this solve
  const bool centred = qualifies && inMiddleNeighbourhood(stats);
  if (!_stored_this_solve || centred || (qualifies && !_middle_centred)) {
    storeCachedPoint(solver, solver._kappa);
    _stored_this_solve = true;
    _middle_centred = _middle_centred || centred;
  }
}

// A solve broke down: linear-solver failure, invalid auxiliaries, a
// non-finite objective, or (not converged and) a final constraint violation
// far above the one it started from.
bool ISPDStartPointPolicy::brokeDown(SparseOptimizerISPD &solver,
                                     bool solverFailed) {
  if (solverFailed || !solver.auxiliariesValid()) {
    return true;
  }
  solver.computeActiveErrors();
  if (!std::isfinite(solver.activeRobustChi2()))
    return true;
  const double factor = _options.recovery.divergenceFactor;
  if (!_keep_rollback || solver._last_converged || factor <= 0.0) {
    return false;
  }
  const double floor = std::max(
      _start_violation, solver.terminationCriteria().epsilon_constraint);
  return !(solver.constraintViolation() <= factor * floor);
}

FinishAction ISPDStartPointPolicy::finish(SparseOptimizerISPD &solver,
                                          const SolveOutcome &outcome) {
  if (outcome.iterations < 0) {
    if (_keep_rollback) {
      restoreRollback(solver);
      _last_restored = true;
    }
    invalidate();
    return FinishAction::Error;
  }

  const bool broken = brokeDown(solver, outcome.solverFailed);
  if (broken && _keep_rollback) {
    // put the starting point back instead of leaving the diverged iterate in
    // the vertices (a middle point cached by this solve is dropped)
    _vertices.discard(VertexEstimateCache::Cached);
    restoreRollback(solver);
    solver._last_converged = false;
    _stored_this_solve = false;
    _last_restored = true;
    if (!_retrying && _started_warm && _options.recovery.retry)
      return FinishAction::RetryCold;
  } else if (_retrying) {
    _last_restored = false; // the cold retry recovered
  }
  if (outcome.solverFailed || broken) {
    // the iterate is not usable: the next solve starts cold
    invalidate();
    return FinishAction::Reject;
  }

  // cache the final iterate (PreviousFinal), or the final one if no middle
  // iterate was cached (PreviousMiddle)
  if (_options.source == ISPDStartSource::PreviousFinal ||
      (_options.source == ISPDStartSource::PreviousMiddle &&
       !_stored_this_solve)) {
    storeCachedPoint(solver, solver._kappa);
    _stored_this_solve = true;
  }
  _warm_point_valid =
      _options.source != ISPDStartSource::Cold && _stored_this_solve;
  _vertices.discard(VertexEstimateCache::Rollback);
  if (!_warm_point_valid)
    _vertices.discard(VertexEstimateCache::Cached);
  return FinishAction::Accept;
}

void ISPDStartPointPolicy::prepareColdRetry(SparseOptimizerISPD &solver) {
  // solve once more from the restored point with a cold start of nu,
  // lambda, s and kappa
  solver.coldStart(true);
  storeRollback(solver);
  _last_source = 0;
  _last_retry = true;
  _retrying = true;
  _middle_centred = false;
}

} // namespace cg2o
