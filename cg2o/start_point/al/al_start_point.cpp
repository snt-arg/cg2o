#include "al_start_point.h"
#include "cg2o/core/sparse_optimizer_al.h"
#include <stdexcept>

namespace cg2o {

void ALStartContext::resetMultipliers() { _start.resetMultipliers(); }

ALStartPoint::ALStartPoint(SparseOptimizerAL &solver)
    : _solver(solver), _method(std::make_unique<al::Cold>()) {}

void ALStartPoint::select(ALStartMethod method) {
  switch (method) {
  case ALStartMethod::Cold:
    use<al::Cold>();
    return;
  case ALStartMethod::ReuseMultipliers:
    use<al::ReuseMultipliers>();
    return;
  default:
    throw std::invalid_argument("AL start method: 0 cold, 1 reuse multipliers");
  }
}

void ALStartPoint::prepare() {
  ALStartContext context(*this);
  _method.active().prepare(context);
}

void ALStartPoint::resetMultipliers() {
  _solver.resetLagrangeMultiplierEq();
  for (auto &edge : _solver._activeEdgesIneq)
    _solver.executeIneqMultiplierFunction(edge, MultiplierOperation::Reset);
}

} // namespace cg2o
