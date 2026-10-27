// Numerical building blocks for ISPD warm starts.
#ifndef CG2O_WARM_START_MATH_H
#define CG2O_WARM_START_MATH_H

#include <algorithm>
#include <cmath>
#include <limits>

namespace cg2o::warm_start {

inline bool positiveFinite(double value) {
  return value > 0.0 && std::isfinite(value);
}

// The quantities used by the eliminated primal-dual system must be finite.
// Positivity of the product alone would incorrectly accept two negatives.
inline bool validPair(double slack, double multiplier, double kappa) {
  return positiveFinite(slack) && positiveFinite(multiplier) &&
         positiveFinite(kappa) && positiveFinite(slack * multiplier) &&
         positiveFinite(multiplier / slack) &&
         positiveFinite((1.0 / kappa) / slack);
}

// Project s*lambda into the requested band about mu. Logs avoid overflow
// in target/(s*lambda). Mode 0 preserves the ratio, 2 preserves slack,
// 3 preserves the multiplier. Mode 1 is the legacy, unit-dependent rule.
// Explicit component floors take precedence if the band is incompatible.
inline bool centerPair(double& slack, double& multiplier, double mu,
                       double lo, double hi, int mode,
                       double slackFloor, double multiplierFloor) {
  if (!positiveFinite(slack) || !positiveFinite(multiplier) ||
      !positiveFinite(mu)) return false;
  const double product = slack * multiplier;
  if (positiveFinite(product) && product >= mu*lo && product <= mu*hi &&
      slack >= slackFloor && multiplier >= multiplierFloor)
    return validPair(slack, multiplier, 1.0/mu);
  const double logS = std::log(slack), logL = std::log(multiplier);
  const double logMu = std::log(mu);
  const double lower = lo > 0.0 ? logMu + std::log(lo)
                               : -std::numeric_limits<double>::infinity();
  const double upper = std::isfinite(hi) ? logMu + std::log(hi)
                                        : std::numeric_limits<double>::infinity();
  const double logP = logS + logL;
  const double smallest = std::numeric_limits<double>::min();
  const double minS = std::log(std::max(slackFloor, smallest));
  const double minL = std::log(std::max(multiplierFloor, smallest));
  if (logP >= lower && logP <= upper && logS >= minS && logL >= minL)
    return true;
  const double target = std::max(std::clamp(logP, lower, upper), minS + minL);
  double newS;
  if (mode == 2 || (mode == 1 && slack >= multiplier)) newS = logS;
  else if (mode == 3 || mode == 1) newS = target - logL;
  else newS = 0.5 * (target + logS - logL);
  newS = std::clamp(newS, minS, std::max(minS, target - minL));
  slack = std::max(std::exp(newS), slackFloor);
  multiplier = std::max(std::exp(target - newS), multiplierFloor);
  return validPair(slack, multiplier, 1.0 / mu);
}

} // namespace cg2o::warm_start
#endif
