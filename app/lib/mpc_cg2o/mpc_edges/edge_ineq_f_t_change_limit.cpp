#include "edge_ineq_f_t_change_limit.h"

namespace cg2o::mpc {

EdgeIneq_f_t_change_limit::EdgeIneq_f_t_change_limit(
    int k, std::shared_ptr<MPCParameters> param = nullptr)
    : _k(k), _param(param) {}

void EdgeIneq_f_t_change_limit::computeIneq() {
  double f_t_k_prev = 0;
  if (this->_k == 0) {
    f_t_k_prev = _param->get_f_t_prev();
  } else {
    const Vertex_f_t *v_f_t_k_prev =
        static_cast<const Vertex_f_t *>(_vertices[0]);
    f_t_k_prev = v_f_t_k_prev->estimate();
  }

  const Vertex_f_t *v_f_t_k = static_cast<const Vertex_f_t *>(_vertices[1]);
  double f_t_k = v_f_t_k->estimate();

  const Vertex_slack_2 *v_slack_2_k =
      static_cast<const Vertex_slack_2 *>(_vertices[2]);
  double slack_2_k = v_slack_2_k->estimate();

  double delta_f_t_max = _param->get_delta_f_t_max();

  _ineq[0] = f_t_k - f_t_k_prev - (delta_f_t_max + slack_2_k);    // g(x) <= 0
  _ineq[1] = -(f_t_k - f_t_k_prev) - (delta_f_t_max + slack_2_k); // g(x) <= 0
}

#ifdef USE_EXACT_JACOBIANS
void EdgeIneq_f_t_change_limit::linearizeOplus() {
  // The derivative of _ineq[0]  with respect to the vertex[0] (f_t_k_prev)  is
  // -1 The derivative of _ineq[1]  with respect to the vertex[0] (f_t_k_prev)
  // is 1 The derivative of _ineq[0]  with respect to the vertex[1] (f_t_k)  is
  // 1 The derivative of _ineq[1]  with respect to the vertex[1] (f_t_k)  is -1
  // The derivative of _ineq[0]  with respect to the vertex[2] (slack_2_k)  is
  // -1 The derivative of _ineq[1]  with respect to the vertex[2] (slack_2_k) is
  // -1 Set the Jacobian's value for the scalar relationship
  auto &J_v0 = std::get<0>(this->_jacobianOplus); // f_t_k_prev
  auto &J_v1 = std::get<1>(this->_jacobianOplus); // f_t_k
  auto &J_v2 = std::get<2>(this->_jacobianOplus); // slack_2_k
  if (_k == 0) {
    J_v0 << 0, // J(e_0,v_0[0])
        0;     // J(e_1,v_0[0])
  } else {
    J_v0 << -1.0, // J(e_0,v_0[0])
        1.0;      // J(e_1,v_0[0])
  }
  J_v1 << 1.0,  // J(e_0,v_1[0])
      -1.0;     // J(e_1,v_1[0])
  J_v2 << -1.0, // J(e_0,v_2[0])
      -1.0;     // J(e_1,v_2[0])
}
#endif

bool EdgeIneq_f_t_change_limit::write(std::ostream &os) const {
  os << "EdgeIneq_f_t_change_limit: ";
  for (int i = 0; i < 2; ++i) {
    os << "ineq[" << i << "] = " << _ineq[i] << " ";
  }
  return os.good();
}

bool EdgeIneq_f_t_change_limit::read(std::istream &is) {
  for (int i = 0; i < 2; ++i) {
    is >> _ineq[i];
  }
  return is.good();
}

} // namespace cg2o::mpc
