#include "bipm_start_point.h"
#include "cg2o/core/sparse_optimizer_bipm.h"
#include <stdexcept>

namespace cg2o {

void BIPMStartContext::coldStart() { _start.coldStart(); }

BIPMStartPoint::BIPMStartPoint(SparseOptimizerBIPM &solver)
    : _solver(solver), _method(std::make_unique<bipm::Cold>()) {}

void BIPMStartPoint::select(BIPMStartMethod method) {
  switch (method) {
  case BIPMStartMethod::Cold:
    use<bipm::Cold>();
    return;
  default:
    throw std::invalid_argument("BIPM start method: 0 cold");
  }
}

void BIPMStartPoint::prepare() {
  BIPMStartContext context(*this);
  _method.active().prepare(context);
}

void BIPMStartPoint::coldStart() {
  // without inequality edges there is no barrier term
  _solver._kappa = _solver._activeEdgesIneq.empty() ? _solver._kappa_final
                                                    : _solver._kappa_initial;
  if (!_options.keepEqualityMultipliers)
    _solver.resetLagrangeMultiplierEq();
}

} // namespace cg2o
