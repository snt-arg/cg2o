#ifndef EDGE_DYNAMICS_D_H_H
#define EDGE_DYNAMICS_D_H_H

#include "cg2o/core/base_fixed_sized_edge_eq.h"
#include "mpc_parameters.h"
#include "mpc_vertices.h"

namespace cg2o::mpc {

class EdgeDynamics_d_h
    : public cg2o::BaseFixedSizedEdgeEq<1, double, Vertex_d_h, Vertex_d_h,
                                       Vertex_v_h, Vertex_v_h> {
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  EdgeDynamics_d_h(int k, std::shared_ptr<MPCParameters> param);
  void computeEq() override;
#ifndef MPC_USE_NUMERICAL_JACOBIAN
  void linearizeOplus() override; // optional
#endif

  bool write(std::ostream &os) const override;
  bool read(std::istream &is) override;

private:
  int _k;
  std::shared_ptr<MPCParameters> _param;
};

} // namespace cg2o::mpc

#endif // EDGE_DYNAMICS_D_H_H
