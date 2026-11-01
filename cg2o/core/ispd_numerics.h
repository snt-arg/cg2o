// Numerical checks on the ISPD auxiliary variables, used by the solver and by
// its start point.
#ifndef CG2O_ISPD_NUMERICS_H
#define CG2O_ISPD_NUMERICS_H

#include <cmath>

namespace cg2o::ispd {

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

} // namespace cg2o::ispd

#endif // CG2O_ISPD_NUMERICS_H
