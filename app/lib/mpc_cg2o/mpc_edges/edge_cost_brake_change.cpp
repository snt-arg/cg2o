#include "edge_cost_brake_change.h"

namespace cg2o::mpc {

EdgeCost_brake_change::EdgeCost_brake_change(int k, std::shared_ptr<MPCParameters> param = nullptr)
    : _k(k), _param(param){}

void EdgeCost_brake_change::computeError() {
    const Vertex_slack_3* v_slack_3_k = static_cast<const Vertex_slack_3*>(_vertices[0]);
    double slack_3_k = v_slack_3_k->estimate();

    _error[0] = slack_3_k;
}

void EdgeCost_brake_change::linearizeOplus() {
    // The derivative of _error[0] (slack_3_k) with respect to the vertex[0] (slack_3_k)  is 1
     // Set the Jacobian's value for the scalar relationship 
    auto& J_v0 = std::get<0>(this->_jacobianOplus);  
    J_v0 << 1.0;   // J(e_0,v_0[0])

}
bool EdgeCost_brake_change::write(std::ostream& os) const {
    os << "EdgeCost_brake_change, ";
    for (int i = 0; i < 1; i++) {
        os << _error[i] << " ";
    }

    return os.good();
}

bool EdgeCost_brake_change::read(std::istream& is) {
    for (int i = 0; i < 1; i++) {
        is >> _error[i];
    }

    return is.good();
    
}
}  // namespace cg2o::mpc
