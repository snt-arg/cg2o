#include "edge_ineq_f_t_min.h"

namespace cg2o::mpc {

EdgeIneq_f_t_min::EdgeIneq_f_t_min(
    int k, std::shared_ptr<MPCParameters> param = nullptr)
    : _k(k), _param(param) {
  if (_param->isScaleInequalities()) {
    _scaling_factor = 10e-5;
  } else {
    _scaling_factor = 1.0;
  }
}

void EdgeIneq_f_t_min::computeIneq() {
  const Vertex_f_t *v_f_t_k = static_cast<const Vertex_f_t *>(_vertices[0]);

  double f_t_k = v_f_t_k->estimate();
  double b3 = _param->get_b3();

  if (_param->isLinearInequalities() || true) { // temp force it true
    _ineq[0] = b3 - f_t_k; // g(x) <= 0 means f_t(k) >= b3
  } else {
    _ineq[0] =
        b3 * b3 * b3 - f_t_k * f_t_k * f_t_k; // g(x) <= 0 means f_t(k) >= b3
    _ineq[0] *= _scaling_factor;
  }
}

#if MPC_USE_NUMERICAL_JACOBIANS == 0
void EdgeIneq_f_t_min::linearizeOplus() {
  // The derivative of _ineq[0]  with respect to the vertex[0] (f_t_k):
  // Set the Jacobian's value for the scalar relationship
  const Vertex_f_t *v_f_t_k = static_cast<const Vertex_f_t *>(_vertices[0]);
  double f_t_k = v_f_t_k->estimate();

  auto &J_v0 = std::get<0>(this->_jacobianOplus);
  if (_param->isLinearInequalities() || true) {
    J_v0 << -1.0; // J(e_0,v_0[0])
  } else {
    J_v0 << -(3.0 * pow(f_t_k, 2)) * _scaling_factor; // J(e_0,v_0[0])
  }
}
#endif
bool EdgeIneq_f_t_min::write(std::ostream &os) const {
  os << "EdgeIneq_f_t_min: ";
  for (int i = 0; i < 1; ++i) {
    os << "ineq[" << i << "] = " << _ineq[i] << " ";
  }
  return os.good();
}

bool EdgeIneq_f_t_min::read(std::istream &is) {
  for (int i = 0; i < 1; ++i) {
    is >> _ineq[i];
  }
  return is.good();
}

} // namespace cg2o::mpc
