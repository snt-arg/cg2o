#include "edge_ineq_f_b_min.h"

namespace cg2o::mpc {

EdgeIneq_f_b_min::EdgeIneq_f_b_min(
    int k, std::shared_ptr<MPCParameters> param = nullptr)
    : _k(k), _param(param) {
  if (_param->isScaleInequalities()) {
    _scaling_factor = 1e-5;
  } else {
    _scaling_factor = 1.0;
  }
}

void EdgeIneq_f_b_min::computeIneq() {
  const Vertex_f_b *v_f_b_k = static_cast<const Vertex_f_b *>(_vertices[0]);

  double f_b_k = v_f_b_k->estimate();

  if (_param->isLinearInequalities() || true) {
    _ineq[0] = 0 - f_b_k; // g(x) <= 0 means f_b(k) >= 0
  } else {
    _ineq[0] = 0 - pow(f_b_k, 3); // g(x) <= 0 means f_b(k) >= 0
    _ineq[0] *= _scaling_factor;
  }
}
#if MPC_USE_NUMERICAL_JACOBIANS == 0
void EdgeIneq_f_b_min::linearizeOplus() {
  const Vertex_f_b *v_f_b_k = static_cast<const Vertex_f_b *>(_vertices[0]);
  double f_b_k = v_f_b_k->estimate();

  // The derivative of _ineq[0]  with respect to the vertex[0] (f_b_k)  is -1
  // Set the Jacobian's value for the scalar relationship
  auto &J_v0 = std::get<0>(this->_jacobianOplus);
  if (_param->isLinearInequalities() || true) {
    J_v0 << -1.0; // J(e_0,v_0[0])
  } else {
    J_v0 << -3.0 * pow(f_b_k, 2) * _scaling_factor; // J(e_0,v_0[0])
  }
}
#endif
bool EdgeIneq_f_b_min::write(std::ostream &os) const {
  os << "EdgeIneq_f_b_min: ";
  for (int i = 0; i < 1; ++i) {
    os << "ineq[" << i << "] = " << _ineq[i] << " ";
  }
  return os.good();
}

bool EdgeIneq_f_b_min::read(std::istream &is) {
  for (int i = 0; i < 1; ++i) {
    is >> _ineq[i];
  }
  return is.good();
}

} // namespace cg2o::mpc
