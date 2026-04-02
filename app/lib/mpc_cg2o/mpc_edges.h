#ifndef MPC_EDGES_H
#define MPC_EDGES_H

// Cost Edges
#include "mpc_edges/edge_cost_brake_change.h"
#include "mpc_edges/edge_cost_brake.h"
#include "mpc_edges/edge_cost_energy.h"
#include "mpc_edges/edge_cost_inner_distance.h"
#include "mpc_edges/edge_cost_traction_change.h"


// Equality Constraints
#include "mpc_edges/edge_dynamics_v_h.h"
#include "mpc_edges/edge_dynamics_d_h.h"

// Inequality Constraints
#include "mpc_edges/edge_ineq_v_h_max.h"
#include "mpc_edges/edge_ineq_v_h_min.h"
#include "mpc_edges/edge_ineq_d_h_max.h"
#include "mpc_edges/edge_ineq_d_h_min.h"
#include "mpc_edges/edge_ineq_f_t_max.h"
#include "mpc_edges/edge_ineq_f_t_min.h"
#include "mpc_edges/edge_ineq_f_b_max.h"
#include "mpc_edges/edge_ineq_f_b_min.h"
#include "mpc_edges/edge_ineq_f_t_change_limit.h"
#include "mpc_edges/edge_ineq_f_b_change_limit.h"



#endif  // MPC_EDGES_H

