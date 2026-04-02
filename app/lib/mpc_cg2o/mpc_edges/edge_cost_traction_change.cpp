#include "edge_cost_traction_change.h"

namespace cg2o::mpc {

EdgeCost_traction_change::EdgeCost_traction_change(
    int k, std::shared_ptr<MPCParameters> param = nullptr)
    : _k(k), _param(param) {}

void EdgeCost_traction_change::computeError() {
  const Vertex_slack_2 *v_slack_2_k =
      static_cast<const Vertex_slack_2 *>(_vertices[0]);
  double slack_2_k = v_slack_2_k->estimate();

  _error[0] = slack_2_k;
}
#ifdef USE_EXACT_JACOBIANS
void EdgeCost_traction_change::linearizeOplus() {
  // The derivative of _error[0] (slack_2_k) with respect to the vertex[0]
  // (slack_2_k)  is 1 Set the Jacobian's value for the scalar relationship
  auto &J_v0 = std::get<0>(this->_jacobianOplus);
  J_v0 << 1.0; // J(e_0,v_0[0])
}
#endif

bool EdgeCost_traction_change::write(std::ostream &os) const {
  os << "EdgeCost_traction_change, ";
  for (int i = 0; i < 1; i++) {
    os << _error[i] << " ";
  }

  return os.good();
}

bool EdgeCost_traction_change::read(std::istream &is) {
  for (int i = 0; i < 1; i++) {
    is >> _error[i];
  }

  return is.good();
}

} // namespace cg2o::mpc
