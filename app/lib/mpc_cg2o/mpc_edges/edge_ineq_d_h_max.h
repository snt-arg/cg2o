#ifndef EDGE_INEQ_D_H_MAX_H
#define EDGE_INEQ_D_H_MAX_H

#include "cg2o/core/base_fixed_sized_edge_ineq.h"
#include "mpc_parameters.h"
#include "mpc_vertices.h"

namespace cg2o::mpc {

class EdgeIneq_d_h_max
    : public cg2o::BaseFixedSizedEdgeIneq<1, double, Vertex_d_h, Vertex_v_h,
                                         Vertex_slack_1> {
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  EdgeIneq_d_h_max(int k, std::shared_ptr<MPCParameters> param);
  void computeIneq() override;
#ifdef USE_EXACT_JACOBIANS
  void linearizeOplus() override; // optional
#endif

  bool write(std::ostream &os) const override;
  bool read(std::istream &is) override;

private:
  [[maybe_unused]] int _k;
  std::shared_ptr<MPCParameters> _param;
};

} // namespace cg2o::mpc

#endif // EDGE_INEQ_D_H_MAX_H
