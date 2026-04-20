#include "edge_ineq_f_b_change_limit.h"

namespace cg2o::mpc {

EdgeIneq_f_b_change_limit::EdgeIneq_f_b_change_limit(
    int k, std::shared_ptr<MPCParameters> param = nullptr)
    : _k(k), _param(param) {}

void EdgeIneq_f_b_change_limit::computeIneq() {
  double f_b_k_prev = 0;

  if (this->_k == 0) {
    f_b_k_prev = _param->get_f_b_prev();
  } else {
    const Vertex_f_b *v_f_b_k_prev =
        static_cast<const Vertex_f_b *>(_vertices[0]);
    f_b_k_prev = v_f_b_k_prev->estimate();
  }

  const Vertex_f_b *v_f_b_k = static_cast<const Vertex_f_b *>(_vertices[1]);
  double f_b_k = v_f_b_k->estimate();

  const Vertex_slack_3 *v_slack_3_k =
      static_cast<const Vertex_slack_3 *>(_vertices[2]);
  double slack_3_k = v_slack_3_k->estimate();

  double delta_f_b_max = _param->get_delta_f_b_max();

  _ineq[0] = f_b_k - f_b_k_prev - (delta_f_b_max + slack_3_k);    // g(x) <= 0
  _ineq[1] = -(f_b_k - f_b_k_prev) - (delta_f_b_max + slack_3_k); // g(x) <= 0
}

#ifndef MPC_USE_NUMERICAL_JACOBIAN
void EdgeIneq_f_b_change_limit::linearizeOplus() {
  // The derivative of _ineq[0]  with respect to the vertex[0] (f_b_k_prev)  is
  // -1 The derivative of _ineq[1]  with respect to the vertex[0] (f_b_k_prev)
  // is 1 The derivative of _ineq[0]  with respect to the vertex[1] (f_b_k)  is
  // 1 The derivative of _ineq[1]  with respect to the vertex[1] (f_b_k)  is -1
  // The derivative of _ineq[0]  with respect to the vertex[2] (slack_3_k)  is
  // -1 The derivative of _ineq[1]  with respect to the vertex[2] (slack_3_k) is
  // -1 Set the Jacobian's value for the scalar relationship
  auto &J_v0 = std::get<0>(this->_jacobianOplus); // f_b_k_prev
  auto &J_v1 = std::get<1>(this->_jacobianOplus); // f_b_k
  auto &J_v2 = std::get<2>(this->_jacobianOplus); // slack_3_k

  if (_k == 0) {
    J_v0 << 0, // J(e_0,v_0[0])
        0;     // J(e_1,v_0[0])
  } else {
    J_v0 << -1, // J(e_0,v_0[0])
        1;      // J(e_1,v_0[0])
  }

  J_v1 << 1.0, // J(e_0,v_1[0])
      -1.0;    // J(e_1,v_1[0])

  J_v2 << -1.0, // J(e_0,v_2[0])
      -1.0;     // J(e_1,v_2[0])
}
#endif

bool EdgeIneq_f_b_change_limit::write(std::ostream &os) const {
  os << "EdgeIneq_f_b_change_limit: ";
  for (int i = 0; i < 2; ++i) {
    os << "ineq[" << i << "] = " << _ineq[i] << " ";
  }
  return os.good();
}

bool EdgeIneq_f_b_change_limit::read(std::istream &is) {
  for (int i = 0; i < 2; ++i) {
    is >> _ineq[i];
  }
  return is.good();
}

} // namespace cg2o::mpc
