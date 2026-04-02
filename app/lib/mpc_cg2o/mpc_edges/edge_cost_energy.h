#ifndef EDGE_COST_ENERGY_H
#define EDGE_COST_ENERGY_H


#include "mpc_parameters.h"
#include "mpc_vertices.h"
#include "g2o/core/base_fixed_sized_edge.h"


namespace cg2o::mpc {

class EdgeCost_energy : public g2o::BaseFixedSizedEdge<2, double, Vertex_v_h, Vertex_f_t> {
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    EdgeCost_energy(int k, std::shared_ptr<MPCParameters> param);
    void computeError() override;
    void linearizeOplus() override; //optional
    
    bool write(std::ostream& os) const override;
    bool read(std::istream& is) override;

    Eigen::Matrix2d get_Q_c();
    Eigen::Vector2d get_b_c();

private:
    [[maybe_unused]] int _k;
    std::shared_ptr<MPCParameters> _param;
   
};

}  // namespace cg2o::mpc

#endif  // EDGE_COST_ENERGY_H
