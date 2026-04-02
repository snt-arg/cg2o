#ifndef EDGE_COST_INNER_DISTANCE_H
#define EDGE_COST_INNER_DISTANCE_H

#include "g2o/core/base_fixed_sized_edge.h"
#include "mpc_parameters.h"
#include "mpc_vertices.h"

namespace cg2o::mpc {

class EdgeCost_inner_distance
    : public g2o::BaseFixedSizedEdge<1, double, Vertex_slack_1> {
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  EdgeCost_inner_distance(int k, std::shared_ptr<MPCParameters> param);
  void computeError() override;
#ifdef USE_EXACT_JACOBIANS
  void linearizeOplus() override;
#endif

  bool write(std::ostream &os) const override;
  bool read(std::istream &is) override;

private:
  [[maybe_unused]] int _k;
  std::shared_ptr<MPCParameters> _param;
};

} // namespace cg2o::mpc

#endif // EDGE_COST_INNER_DISTANCE_H
