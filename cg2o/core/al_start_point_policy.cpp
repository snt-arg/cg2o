#include "al_start_point_policy.h"
#include "sparse_optimizer_al.h"

namespace cg2o {

void ALStartPointPolicy::prepare(SparseOptimizerAL &solver) {
  if (_options.source == ALStartSource::Cold) {
    // initial multipliers and penalties of the equality and inequality edges
    solver._reset_lagrange_multipliers = true;
    solver.resetLagrangeMultiplierEq();
    for (auto &edge : solver._activeEdgesIneq) {
      solver.executeIneqMultiplierUpdate(edge);
    }
  }
  solver._reset_lagrange_multipliers = false;
}

} // namespace cg2o
