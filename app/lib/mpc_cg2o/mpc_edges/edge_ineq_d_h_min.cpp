#include "edge_ineq_d_h_min.h"
#include <cmath>

namespace cg2o::mpc {

EdgeIneq_d_h_min::EdgeIneq_d_h_min(
    int k, std::shared_ptr<MPCParameters> param = nullptr)
    : _k(k), _param(param) {}

void EdgeIneq_d_h_min::computeIneq() {
  const Vertex_d_h *v_d_h_kp1 = static_cast<const Vertex_d_h *>(_vertices[0]);
  const Vertex_v_h *v_v_h_kp1 = static_cast<const Vertex_v_h *>(_vertices[1]);

  double d_h_kp1 = v_d_h_kp1->estimate();
  double v_h_kp1 = v_v_h_kp1->estimate();

  double h_min = _param->get_h_min();
  double d_min = _param->get_d_min();

  if (_param->isLinearInequalities() || true) {
    _ineq[0] = (h_min * v_h_kp1 + d_min) - d_h_kp1;
  } else {
    _ineq[0] = pow(h_min * v_h_kp1 + d_min, 3) - pow(d_h_kp1, 3);
  }
}
#if MPC_USE_NUMERICAL_JACOBIANS == 0
void EdgeIneq_d_h_min::linearizeOplus() {
  // The derivative of _ineq[0]  with respect to the vertex[0] (d_h_kp1)  is -1
  // The derivative of _ineq[0]  with respect to the vertex[1] (v_h_kp1)  is
  // h_min Set the Jacobian's value for the scalar relationship
  auto &J_v0 = std::get<0>(this->_jacobianOplus); // d_h_kp1
  auto &J_v1 = std::get<1>(this->_jacobianOplus); // v_h_kp1

  if (_param->isLinearInequalities() || true) {
    double h_min = _param->get_h_min();

    J_v0 << -1.0;  // J(e_0,v_0[0])
    J_v1 << h_min; // J(e_0,v_1[0])
    return;
  } else {
    const Vertex_d_h *v_d_h_kp1 = static_cast<const Vertex_d_h *>(_vertices[0]);
    const Vertex_v_h *v_v_h_kp1 = static_cast<const Vertex_v_h *>(_vertices[1]);

    double d_h_kp1 = v_d_h_kp1->estimate();
    double v_h_kp1 = v_v_h_kp1->estimate();
    double h_min = _param->get_h_min();
    double d_min = _param->get_d_min();

    J_v0 << -3 * pow(d_h_kp1, 2);                        // J(e_0,v_0[0])
    J_v1 << 3 * pow(h_min * v_h_kp1 + d_min, 2) * h_min; // J(e_0,v_1[0])
    return;
  }
}

#endif

bool EdgeIneq_d_h_min::write(std::ostream &os) const {
  os << "EdgeIneq_d_h_min, ";
  for (int i = 0; i < 1; i++) {
    os << _ineq[i] << " ";
  }

  return os.good();
}

bool EdgeIneq_d_h_min::read(std::istream &is) {
  for (int i = 0; i < 1; i++) {
    is >> _ineq[i];
  }

  return is.good();
}

} // namespace cg2o::mpc
