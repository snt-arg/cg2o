#ifndef EDGE_INEQ_F_T_CHANGE_LIMIT_H
#define EDGE_INEQ_F_T_CHANGE_LIMIT_H

#include "cg2o/core/base_fixed_sized_edge_ineq.h"
#include "mpc_parameters.h"
#include "mpc_vertices.h"

namespace cg2o::mpc {

class EdgeIneq_f_t_change_limit
    : public cg2o::BaseFixedSizedEdgeIneq<2, double, Vertex_f_t, Vertex_f_t,
                                         Vertex_slack_2> {
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  EdgeIneq_f_t_change_limit(int k, std::shared_ptr<MPCParameters> param);
  void computeIneq() override;
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

#endif // EDGE_INEQ_F_T_CHANGE_LIMIT_1_H
