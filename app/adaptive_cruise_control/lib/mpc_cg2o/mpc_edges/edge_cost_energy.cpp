#include "edge_cost_energy.h"
#include <Eigen/Dense>

namespace cg2o::mpc {

EdgeCost_energy::EdgeCost_energy(int k, std::shared_ptr<MPCParameters> param = nullptr)
    : _k(k), _param(param){}

void EdgeCost_energy::computeError() {
    const Vertex_v_h* v_v_h_k = static_cast<const Vertex_v_h*>(_vertices[0]);
    const Vertex_f_t* v_f_t_k = static_cast<const Vertex_f_t*>(_vertices[1]);

    double v_h_k = v_v_h_k->estimate();
    double f_t_k = v_f_t_k->estimate();

    Eigen::Vector2d temp_e(v_h_k * v_h_k, f_t_k);

    Eigen::Vector2d b_c = get_b_c();

    Eigen::Matrix2d Q_c = get_Q_c();
  
    Eigen::Matrix2d Q_c_inv = Q_c.inverse();

    
    _error = temp_e + 0.5 * Q_c_inv * b_c ;
}


    bool EdgeCost_energy::write(std::ostream& os) const {
    os << "EdgeCost_brake_change, ";
    for (int i = 0; i < 1; i++) {
        os << _error[i] << " ";
    }

    return os.good();

}

Eigen::Matrix2d EdgeCost_energy::get_Q_c() {

    std::vector<double> eff_map = _param->get_eff_map();  
    //double p01 = eff_map[1];
    double p02 = eff_map[2];
    // double p10 = eff_map[3];
    double p11 = eff_map[4];
    double p20 = eff_map[5];   

 
    Eigen::Matrix2d Q_c;
    Q_c << 2 * p20, p11,
               p11 , 2 * p02; // this 2 Q_c therefor we multiply by 0.5 
    Q_c = .5 * Q_c;
    return Q_c;
}

void EdgeCost_energy::linearizeOplus() {
    auto& J_v0 = std::get<0>(this->_jacobianOplus);  
    auto& J_v1 = std::get<1>(this->_jacobianOplus); 
   
    const Vertex_v_h* v_v_h_k = static_cast<const Vertex_v_h*>(_vertices[0]);    
    double v_h_k = v_v_h_k->estimate();
 
    J_v0 << 2 * v_h_k,       // J(e0,v0[0]), ..., J(e0,v0[v0.size() - 1])
            0.0;       // J(e1,v0[0]), ..., J(e1,v0[v0.size() - 1])
    J_v1 << 0.0,       // J(e0,v1[0]), ..., J(e0,v0[v1.size() - 1])
            1.0;       // J(e1,v1[0]), ..., J(e1,v0[v1.size() - 1])
}

Eigen::Vector2d EdgeCost_energy::get_b_c() {
    std::vector<double> eff_map = _param->get_eff_map();  
    double p01 = eff_map[1];
    double p10 = eff_map[3];


    Eigen::Vector2d b_c(p10, p01);
    return b_c;
}

bool EdgeCost_energy::read(std::istream& is) {
    for (int i = 0; i < 1; i++) {
        is >> _error[i];
    }

    return is.good();
    
}
}  // namespace cg2o::mpc
