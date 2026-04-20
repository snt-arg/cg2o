#ifndef EDGE_COST_TRACTION_CHANGE_H
#define EDGE_COST_TRACTION_CHANGE_H

#include "g2o/core/base_fixed_sized_edge.h"
#include "mpc_parameters.h"
#include "mpc_vertices.h"

namespace cg2o::mpc {

class EdgeCost_traction_change
    : public g2o::BaseFixedSizedEdge<1, double, Vertex_slack_2> {
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  EdgeCost_traction_change(int k, std::shared_ptr<MPCParameters> param);
  void computeError() override;
#ifndef MPC_USE_NUMERICAL_JACOBIAN
  void linearizeOplus() override; // optional
#endif

  bool write(std::ostream &os) const override;
  bool read(std::istream &is) override;

private:
  [[maybe_unused]] int _k;
  std::shared_ptr<MPCParameters> _param;
};

} // namespace cg2o::mpc

#endif // EDGE_COST_TRACTION_CHANGE_H
