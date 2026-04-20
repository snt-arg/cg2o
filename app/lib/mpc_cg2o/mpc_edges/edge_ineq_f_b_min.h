#ifndef EDGE_INEQ_F_B_MIN_H
#define EDGE_INEQ_F_B_MIN_H

#include "cg2o/core/base_fixed_sized_edge_ineq.h"
#include "mpc_parameters.h"
#include "mpc_vertices.h"

namespace cg2o::mpc {

class EdgeIneq_f_b_min
    : public cg2o::BaseFixedSizedEdgeIneq<1, double, Vertex_f_b> {
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  EdgeIneq_f_b_min(int k, std::shared_ptr<MPCParameters> param);
  void computeIneq() override;
#ifndef MPC_USE_NUMERICAL_JACOBIAN
  void linearizeOplus() override; // optional
#endif
  bool write(std::ostream &os) const override;
  bool read(std::istream &is) override;

private:
  [[maybe_unused]] int _k;
  std::shared_ptr<MPCParameters> _param;
  double _scaling_factor;
};

} // namespace cg2o::mpc

#endif // EDGE_INEQ_F_B_MIN_H
