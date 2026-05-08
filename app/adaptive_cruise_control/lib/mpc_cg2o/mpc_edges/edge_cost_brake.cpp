#include "edge_cost_brake.h"

namespace cg2o::mpc {

EdgeCost_brake::EdgeCost_brake(int k, std::shared_ptr<MPCParameters> param = nullptr)
    : _k(k), _param(param){}

void EdgeCost_brake::computeError() {
    const Vertex_f_b* v_f_b_k = static_cast<const Vertex_f_b*>(_vertices[0]);

    double f_b_k = v_f_b_k->estimate();
 
    _error[0] = f_b_k;
}

void EdgeCost_brake::linearizeOplus() {
    // The derivative of _error[0] (f_b_k) with respect to the vertex[0] (f_b_k)  is 1
     // Set the Jacobian's value for the scalar relationship 
    auto& J_v0 = std::get<0>(this->_jacobianOplus);  
    J_v0 << 1.0;   // J(e_0,v_0[0])

}


bool EdgeCost_brake::write(std::ostream& os) const {
    os << "EdgeCost_brake_change, ";
    for (int i = 0; i < 1; i++) {
        os << _error[i] << " ";
    }

    return os.good();

}

bool EdgeCost_brake::read(std::istream& is) {
    for (int i = 0; i < 1; i++) {
        is >> _error[i];
    }

    return is.good();
    
}
}  // namespace cg2o::mpc
