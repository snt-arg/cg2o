#include "ispd_strategies.h"
#include "cg2o/core/ispd_numerics.h"
#include <cmath>
#include <stdexcept>

namespace cg2o::ispd {

namespace {
// the stored point if it can be used, a cold start otherwise
std::optional<double> storedOrCold(ISPDStartContext &start) {
  if (start.storedPointUsable())
    return start.restoreStoredPoint();
  start.coldStart();
  return std::nullopt;
}
} // namespace

std::optional<double> Cold::load(ISPDStartContext &start) {
  start.coldStart();
  return std::nullopt;
}

bool PreviousFinal::available(const ISPDStartContext &start) const {
  return start.storedPointUsable();
}

std::optional<double> PreviousFinal::load(ISPDStartContext &start) {
  return storedOrCold(start);
}

void PreviousFinal::retain(ISPDStartContext &start) {
  start.storeCurrentPoint();
}

void PreviousMiddle::validate() const {
  const auto &m = settings;
  const int rule = static_cast<int>(m.rule);
  if (rule < 0 || rule > 3 || !positiveFinite(m.kappa) || m.iteration < 0 ||
      !positiveFinite(m.gap) ||
      !(m.neighbourhoodLo >= 0.0 && m.neighbourhoodLo <= 1.0 &&
        m.neighbourhoodHi >= 1.0))
    throw std::invalid_argument("invalid PreviousMiddle settings");
}

bool PreviousMiddle::available(const ISPDStartContext &start) const {
  return start.storedPointUsable();
}

std::optional<double> PreviousMiddle::load(ISPDStartContext &start) {
  return storedOrCold(start);
}

IterateObservation PreviousMiddle::observation() const {
  // the relative gap needs the objective of each iterate
  return settings.rule == MiddlePointRule::RelativeGap
             ? IterateObservation::ComplementarityAndCost
             : IterateObservation::Complementarity;
}

bool PreviousMiddle::inNeighbourhood(const ISPDPointStats &stats) const {
  return stats.rows == 0 || (stats.minRatio >= settings.neighbourhoodLo &&
                             stats.maxRatio <= settings.neighbourhoodHi);
}

void PreviousMiddle::observe(ISPDStartContext &start,
                             const ISPDPointStats &stats) {
  // keep the last iterate that qualifies as the middle point (the first
  // iterate if none does)
  if (!(std::isfinite(stats.mu) && stats.mu > 0.0 && start.pairsValid()))
    return;
  const auto &m = settings;
  bool qualifies = false;
  switch (m.rule) {
  case MiddlePointRule::BarrierKappa:
    qualifies = start.kappa() <= m.kappa;
    break;
  case MiddlePointRule::DualityGapKappa:
    qualifies = (stats.mu > 0.0 ? 1.0 / stats.mu : 0.0) <= m.kappa;
    break;
  case MiddlePointRule::Iteration:
    qualifies = stats.iteration < m.iteration;
    break;
  case MiddlePointRule::RelativeGap:
    qualifies = stats.relativeGap >= m.gap;
    break;
  default:
    throw std::runtime_error("[Error] Unknown middle point rule");
  }
  // the neighbourhood is a preference: a qualifying iterate inside it
  // replaces any stored one, a qualifying iterate outside it is stored
  // only while no centred one has been stored in this attempt
  const bool centred = qualifies && inNeighbourhood(stats);
  if (!start.storedInThisSolve() || centred || (qualifies && !_centred)) {
    start.storeCurrentPoint();
    _centred = _centred || centred;
  }
}

void PreviousMiddle::retain(ISPDStartContext &start) {
  // no iterate was stored during the solve: the final one
  if (!start.storedInThisSolve())
    start.storeCurrentPoint();
}

std::optional<double> CurrentGraph::load(ISPDStartContext &start) {
  return start.takeGraphPoint();
}

} // namespace cg2o::ispd
