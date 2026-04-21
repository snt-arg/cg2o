#include "edge_cost_inner_distance.h"

namespace cg2o::mpc {

EdgeCost_inner_distance::EdgeCost_inner_distance(
    int k, std::shared_ptr<MPCParameters> param = nullptr)
    : _k(k), _param(param) {}

void EdgeCost_inner_distance::computeError() {
  const Vertex_slack_1 *v_slack_1_kp1 =
      static_cast<const Vertex_slack_1 *>(_vertices[0]);
  double slack_1_kp1 = v_slack_1_kp1->estimate();

  _error[0] = slack_1_kp1;
}

#if MPC_USE_NUMERICAL_JACOBIANS == 0
void EdgeCost_inner_distance::linearizeOplus() {
  // The derivative of _error[0] (slack_1_kp1) with respect to the vertex[0]
  // (slack_1_kp1)  is 1 Set the Jacobian's value for the scalar relationship
  auto &J_v0 = std::get<0>(this->_jacobianOplus);
  J_v0 << 1.0; // J(e_0,v_0[0])
}
#endif

bool EdgeCost_inner_distance::write(std::ostream &os) const {
  os << "EdgeCost_inner_distance, ";
  for (int i = 0; i < 1; i++) {
    os << _error[i] << " ";
  }

  return os.good();
}

bool EdgeCost_inner_distance::read(std::istream &is) {
  for (int i = 0; i < 1; i++) {
    is >> _error[i];
  }

  return is.good();
}

} // namespace cg2o::mpc
