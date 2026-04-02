#ifndef MPC_FORMULATION_H
#define MPC_FORMULATION_H

#include "cg2o/core/base_fixed_sized_edge_eq.h"
#include "cg2o/core/base_fixed_sized_edge_ineq.h"
#include "cg2o/core/sparse_optimizer_bipm.h"
#include "cg2o/core/sparse_optimizer_for_constraints.h"
#include "g2o/core/base_fixed_sized_edge.h"
#include "mpc_edges.h"
#include "mpc_vertices.h"
#include "sparse_optimizer.h"
#include <memory>
#include <optimizable_graph.h>
#include <valarray>
#include <vector>

namespace cg2o::mpc {

template <typename Optimizer> class MPCFormulation {
public:
  MPCFormulation(int N, std::shared_ptr<Optimizer> optimizer,
                 std::shared_ptr<MPCParameters> param);

  // Add vertices
  void addVertex(int offset, int k);
  void addVertices(int k);
  void setVertexScalarEstimate(int id, double value);
  double getVertexScalarEstimate(int id);
  std::vector<double> getResults();
  double getForceInput(int k);
  double getVelocityDesired();
  double getAccelerationDesired(double v_h, double road_grad = 0.0);

  // Add edges
  void addEdgeCost_brake_change(int k);
  void addEdgeCost_brake(int k);
  void addEdgeCost_energy(int k);
  void addEdgeCost_inner_distance(int k);
  void addEdgeCost_traction_change(int k);
  void addEdgeDynamics_v_h(int k);
  void addEdgeDynamics_d_h(int k);
  void addEdgeIneq_v_h_max(int k);
  void addEdgeIneq_v_h_min(int k);
  void addEdgeIneq_d_h_max(int k);
  void addEdgeIneq_d_h_min(int k);
  void addEdgeIneq_f_t_max(int k);
  void addEdgeIneq_f_t_min(int k);
  void addEdgeIneq_f_b_max(int k);
  void addEdgeIneq_f_b_min(int k);
  void addEdgeIneq_f_t_change_limit(int k);
  void addEdgeIneq_f_b_change_limit(int k);
  void addEdges(int k);

  void setN(int N);
  int getN() const;
  void setVLagrangianInitial(double vInitial);
  double getVLagrangianInitial();

  // MPC setup
  void setupMPC();

  double computeCost();

  // Set the initial guess for the optimization
  bool setInitialGuess(bool reset_lagrange_multipliers = true);
  bool setFixedInitialGuess(double value);
  void computeOffsets();

  // For closed loop simulation
  double simulatorVelocity(double u_car);
  void sensor(const double &v_p, std::deque<double> &a_p_memory); 
  std::vector<double> predictor(double v_p,
                                const std::deque<double> &a_p_memory);

protected:
  int _N;
  std::shared_ptr<Optimizer> _optimizer;
  std::shared_ptr<MPCParameters> _param;
  std::vector<double> _results; // v_d[0,N], d_h[0,N], f_t[0,N-1], f_b[0,N-1],
                                // slack_1[1,N], slack_2[0,N-1], slack_3[0,N-1]
  double _init_lagrange_multiplier =
      0.0;                    // Initial value for the Lagrangian vertex
  double _replace_nan = 10.0; // Value to replace NaN estimates

  // Offsets for easy ID calculation
  int _offset_d_h;
  int _offset_slack_1;
  int _offset_f_t;
  int _offset_slack_2;
  int _offset_v_h;
  int _offset_f_b;
  int _offset_slack_3;
  int _offset_eqDynamics_v_h; // for the Lagrangian variables of the equality
                              // constraints
  int _offset_eqDynamics_d_h; // for the Lagrangian variables of the equality
                              // constraints

  // Vertices
  std::vector<std::shared_ptr<g2o::OptimizableGraph::Vertex>> _vertices;
  std::vector<std::shared_ptr<g2o::OptimizableGraph::Edge>> _edges;
  std::vector<std::shared_ptr<g2o::OptimizableGraph::Edge>> _edges_cost;
};

} // namespace cg2o::mpc

// Include the implementation
#include "mpc_formulation.hpp"

#endif // MPC_FORMULATION_H
