#include "edge_ineq_d_h_max.h"
#include "mpc_vertices.h"

#include <cmath>

namespace cg2o::mpc {

EdgeIneq_d_h_max::EdgeIneq_d_h_max(
    int k, std::shared_ptr<MPCParameters> param = nullptr)
    : _k(k), _param(param) {}

void EdgeIneq_d_h_max::computeIneq() {
  const Vertex_d_h *v_d_h_kp1 = static_cast<const Vertex_d_h *>(_vertices[0]);
  const Vertex_v_h *v_v_h_kp1 = static_cast<const Vertex_v_h *>(_vertices[1]);
  const Vertex_slack_1 *v_slack_1 =
      static_cast<const Vertex_slack_1 *>(_vertices[2]);

  double d_h_kp1 = v_d_h_kp1->estimate();
  double v_h_kp1 = v_v_h_kp1->estimate();
  double slack_1_kp1 = v_slack_1->estimate();

  double h_soft = _param->get_h_soft();

  if (_param->isLinearInequalities() || true) {
    _ineq[0] = d_h_kp1 - (h_soft * v_h_kp1 + slack_1_kp1);
  } else {
    _ineq[0] = pow(d_h_kp1, 3) - pow(h_soft * v_h_kp1 + slack_1_kp1, 3);
  }
}
#ifdef USE_EXACT_JACOBIANS
void EdgeIneq_d_h_max::linearizeOplus() {
  auto &J_v0 = std::get<0>(this->_jacobianOplus); // d_h_kp1
  auto &J_v1 = std::get<1>(this->_jacobianOplus); // v_h_kp1
  auto &J_v2 = std::get<2>(this->_jacobianOplus); // slack_1_kp1

  if (_param->isLinearInequalities() || true) {
    J_v0 << 1.0;                   // J(e_0,v_0[0])
    J_v1 << -_param->get_h_soft(); // J(e_0,v_1[0])
    J_v2 << -1.0;                  // J(e_0,v_2[0])
    return;
  } else {
    const Vertex_d_h *v_d_h_kp1 = static_cast<const Vertex_d_h *>(_vertices[0]);
    const Vertex_v_h *v_v_h_kp1 = static_cast<const Vertex_v_h *>(_vertices[1]);
    const Vertex_slack_1 *v_slack_1 =
        static_cast<const Vertex_slack_1 *>(_vertices[2]);

    double d_h_kp1 = v_d_h_kp1->estimate();
    double v_h_kp1 = v_v_h_kp1->estimate();
    double slack_1_kp1 = v_slack_1->estimate();
    double h_soft = _param->get_h_soft();
    J_v0 << 3 * pow(d_h_kp1, 2); // J(e_0,v_0[0])
    J_v1 << -3 * pow(h_soft * v_h_kp1 + slack_1_kp1, 2) *
                h_soft;                                  // J(e_0,v_1[0])
    J_v2 << -3 * pow(h_soft * v_h_kp1 + slack_1_kp1, 2); // J(e_0,v_2[0])
    return;
  }
}
#endif

bool EdgeIneq_d_h_max::write(std::ostream &os) const {
  os << "EdgeIneq_d_h_max, ";
  for (int i = 0; i < 1; i++) {
    os << _ineq[i] << " ";
  }

  return os.good();
}

bool EdgeIneq_d_h_max::read(std::istream &is) {
  for (int i = 0; i < 1; i++) {
    is >> _ineq[i];
  }

  return is.good();
}

} // namespace cg2o::mpc
