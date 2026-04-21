#include "edge_ineq_v_h_min.h"
#include <cmath>

namespace cg2o::mpc {

EdgeIneq_v_h_min::EdgeIneq_v_h_min(
    int k, std::shared_ptr<MPCParameters> param = nullptr)
    : _k(k), _param(param) {}

void EdgeIneq_v_h_min::computeIneq() {
  const Vertex_v_h *v_v_kp1 = static_cast<const Vertex_v_h *>(_vertices[0]);

  double v_h_kp1 = v_v_kp1->estimate();
  double v_min_kp1 = _param->get_v_min_prediction(_k + 1);

  if (_param->isLinearInequalities()) {
    _ineq[0] = v_min_kp1 - v_h_kp1; // g(x) <= 0 means v_h(k+1) >= v_min(k+1)
  } else {
    _ineq[0] = pow(v_min_kp1, 3) -
               pow(v_h_kp1, 3); // g(x) <= 0 means v_h(k+1)^3 >= v_min(k+1)^3
  }
}

#if MPC_USE_NUMERICAL_JACOBIANS == 0
void EdgeIneq_v_h_min::linearizeOplus() {
  // The derivative of _ineq[0]  with respect to the vertex[0] (v_h_kp1)  is -1
  // Set the Jacobian's value for the scalar relationship
  const Vertex_v_h *v_v_kp1 = static_cast<const Vertex_v_h *>(_vertices[0]);
  double v_h_kp1 = v_v_kp1->estimate();

  auto &J_v0 = std::get<0>(this->_jacobianOplus);

  if (_param->isLinearInequalities()) {
    J_v0 << -1.0; // J(e_0,v_0[0])
  } else {
    J_v0 << -3.0 * pow(v_h_kp1, 2); // J(e_0,v_0[0])
  }
}
#endif

bool EdgeIneq_v_h_min::write(std::ostream &os) const {
  os << "EdgeIneq_v_h_min: ";
  for (int i = 0; i < 1; ++i) {
    os << "ineq[" << i << "] = " << _ineq[i] << " ";
  }
  return os.good();
}

bool EdgeIneq_v_h_min::read(std::istream &is) {
  for (int i = 0; i < 1; ++i) {
    is >> _ineq[i];
  }
  return is.good();
}

} // namespace cg2o::mpc
