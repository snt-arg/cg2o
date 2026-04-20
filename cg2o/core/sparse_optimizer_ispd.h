/*
Copyright (c) 2023, University of Luxembourg
All rights reserved.

Redistributions and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its
   contributors may be used to endorse or promote products derived from
   this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS 'AS IS'
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
*/

#ifndef G2O_GRAPH_OPTIMIZER_PDIS_H
#define G2O_GRAPH_OPTIMIZER_PDIS_H

#include "g2o/core/g2o_core_api.h"
#include "sparse_optimizer_for_constraints.h"
#include <algorithm>
#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <optimizable_graph.h>
#include <vector>

namespace cg2o {

// forwad declaration
template <int D, typename E, typename... VertexTypes>
class BaseFixedSizedEdgeIneq;

template <int D, typename E, typename... VertexTypes>
class BaseFixedSizedEdgeEq;

// Specialized optimizer inheriting from SparseOptimizerForConstraints
class G2O_CORE_API SparseOptimizerISPD
    : public SparseOptimizerForConstraints<SparseOptimizerISPD> {
  /**
   * @class SparseOptimizerISPD
   * @brief Infeasible-Start Primal–Dual Interior-Point solver for constrained
   * factor graph optimization.
   *
   * This class extends `SparseOptimizer` to support constrained optimization
   * problems within the factor graph framework using an infeasible-start
   * primal–dual interior-point method (ISPD-IPM).
   *
   * Unlike conventional approaches based on penalty or augmented Lagrangian
   * methods, this implementation directly incorporates inequality constraints
   * into the factor graph through a principled formulation derived from the
   * Karush–Kuhn–Tucker (KKT) conditions. :contentReference[oaicite:0]{index=0}
   *
   * The solver enables joint optimization of:
   * - primal variables (graph states),
   * - Lagrange multipliers,
   * - slack variables,
   * - and barrier parameters,
   * within a unified optimization loop, eliminating the need for nested
   * outer/inner iterations and strictly feasible initialization.
   * :contentReference[oaicite:1]{index=1}
   *
   * --------------------------------------------------------------------------
   * Key Features
   * --------------------------------------------------------------------------
   *
   * - **Infeasible-start formulation**:
   *   The solver does not require the initial solution to satisfy inequality
   *   constraints. Feasibility is achieved progressively using slack variables.
   *
   * - **KKT-based inequality factors**:
   *   Each inequality constraint is reformulated as a dedicated factor node
   *   derived from the primal–dual KKT system, allowing seamless integration
   *   into the global sparse system alongside standard cost and equality
   * factors.
   *
   * - **Unified primal–dual updates**:
   *   Primal variables, Lagrange multipliers, slack variables, and barrier
   *   parameters are updated simultaneously within a single optimization loop.
   *
   * - **Preservation of sparsity**:
   *   The formulation maintains compatibility with the factor graph structure,
   *   enabling efficient sparse linear algebra and scalability to large
   * problems.
   *
   * - **Flexible strategy configuration**:
   *   Multiple strategies are exposed for:
   *     - barrier initialization and update,
   *     - slack and multiplier initialization,
   *     - step-size selection,
   *     - inequality prediction,
   *     - auxiliary variable updates.
   *
   *
   * --------------------------------------------------------------------------
   * Optimization Framework
   * --------------------------------------------------------------------------
   *
   * The constrained optimization problem is formulated as:
   *
   *     minimize   Σ ||e_j(X)||²
   *     subject to h_j(X) = 0
   *                g_j(X) ≤ 0
   *
   * Using the ISPD-IPM formulation, inequality constraints are transformed via
   * slack variables:
   *
   *     g_j(X) + s_j = 0,   s_j ≥ 0
   *
   * and the complementarity condition is relaxed as:
   *
   *     s_j ∘ λ_j = 1 / κ
   *
   * where:
   *   - λ_j are Lagrange multipliers,
   *   - s_j are slack variables,
   *   - κ is the barrier parameter.
   *
   * These conditions are linearized and reduced to a factor-graph-compatible
   * system, enabling inequality constraints to contribute directly to the
   * global system matrix.
   *
   *
   * --------------------------------------------------------------------------
   * Barrier Parameter (κ)
   * --------------------------------------------------------------------------
   *
   * The barrier parameter controls the trade-off between objective minimization
   * and constraint enforcement.
   *
   * - Initialized either as a constant or from the surrogate duality gap.
   * - Updated adaptively based on the average complementarity measure.
   * - Gradually increased to tighten constraint satisfaction.
   *
   * The update rule follows:
   *
   *     κ ← max(κ_lower_bound, ν / η̂)
   *
   * where η̂ is the surrogate duality gap and ν is a scaling factor.
   *
   *
   * --------------------------------------------------------------------------
   * Step-Size and Backtracking
   * --------------------------------------------------------------------------
   *
   * Two independent step-size mechanisms are used:
   *
   * - **Inequality feasibility step**:
   *   Ensures that already satisfied constraints remain feasible.
   *
   * - **Auxiliary variable step**:
   *   Ensures positivity of slack variables and multipliers.
   *
   * The final step size is selected using configurable strategies:
   *   - conservative (min),
   *   - averaged,
   *   - inequality-driven,
   *   - auxiliary-driven,
   *   - or full-step.
   *
   *
   * --------------------------------------------------------------------------
   * Auxiliary Variable Handling
   * --------------------------------------------------------------------------
   *
   * Slack variables and Lagrange multipliers are maintained strictly positive
   * using:
   *
   * - backtracking line search, and/or
   * - safeguard strategies:
   *     - scaling strategy,
   *     - correction-to-threshold strategy.
   *
   *
   * --------------------------------------------------------------------------
   * Implementation Notes
   * --------------------------------------------------------------------------
   *
   * - Inequality constraints are implemented as dedicated factor types whose
   *   contributions are accumulated into the global system similarly to
   *   standard cost and equality factors.
   *
   * - Auxiliary variables (slack and multipliers) are stored locally and are
   *   not part of the graph state vector, but they influence the optimization
   *   through the factor construction.
   *
   * - The solver overrides the standard g2o optimization loop to integrate
   *   primal–dual updates, barrier adaptation, and constraint handling.
   *
   *
   * --------------------------------------------------------------------------
   * Applicability
   * --------------------------------------------------------------------------
   *
   * The solver is applicable to:
   * - constrained SLAM and estimation problems,
   * - trajectory optimization,
   * - model predictive control (MPC),
   * - real-time robotic control systems.
   *
   * When no constraints are present, the solver behaves equivalently to the
   * standard g2o optimizer.
   *
   *
   * @note
   * Proper initialization of the graph and variables is required before calling
   * `initializeOptimization()` or `optimize()`.
   *
   */

public:
  SparseOptimizerISPD();  // Default constructor
  ~SparseOptimizerISPD(); // Virtual destructor

  // Function to construct the quadratic form implementation for Augmented
  // Lagrangian Eqaulity Algorithm
  template <int D, typename E, typename... VertexTypes>
  void
  constructQuadraticFormEq(BaseFixedSizedEdgeEq<D, E, VertexTypes...> &edge);

  template <int D, typename E, typename... VertexTypes>
  bool addEdgeEqImpl(BaseFixedSizedEdgeEq<D, E, VertexTypes...> *e);
  /**
   * add the inequlaity edge to the optimizer and the set of equality edges
   * @param e: the edge to be added
   * @returns false if somethings goes wrong
   */

  void resetLagrangeMultiplierEq() override;

  // Function to construct the quadratic form implementation
  template <int D, typename E, typename... VertexTypes>
  void constructQuadraticFormIneq(
      BaseFixedSizedEdgeIneq<D, E, VertexTypes...> &edge);

  template <int D, typename E, typename... VertexTypes>
  bool addEdgeIneqImpl(BaseFixedSizedEdgeIneq<D, E, VertexTypes...> *e);
  /**
   * add the inequlaity edge to the optimizer and the set of inequality edges
   * @param e: the edge to be added
   * @returns false if somethings goes wrong
   */

  template <int D, typename E, typename... VertexTypes, typename EdgeType>
  double getDualityGap(EdgeType *edge);

  template <int D, typename E, typename... VertexTypes, typename EdgeType>
  void edgeProcessing(EdgeType *edge, int controller);

  template <int D, typename E, typename... VertexTypes>
  void
  initializeSlackVariable(BaseFixedSizedEdgeIneq<D, E, VertexTypes...> *edge);

  void executeEdgeProcessing(void *edgePtr, int controller);
  virtual void update(const double *update) override;
  virtual int optimize(int iterations, bool online = false) override;
  void computeKappa(bool input = false);
  double constraintsBacktracking(const double *update);

  // solver parameter

public:
  //

  void setStepSizeStrategy(int strategy);
  void setIneqBacktrackingStepMin(double value);
  void setAuxBacktrackingStepMin(double value);

  void setInitKappaStrategy(int strategy);
  void setKappaInitial(double value);

  void setUpdateKappaStrategy(int strategy);
  void setKappaFinal(double value);
  void setKappaUpdateFactor(double value);
  void setTau(double value);
  void setLimitKappaFinal(bool value);

  void setInitSlackStrategy(int strategy);
  void setSlackVariableInitialIneq(double value);

  void setInitLagrangeStrategy(int strategy);
  void setLagrangeMultiplierInitialIneq(double value);

  void setIneqPredictionStrategy(int strategy);
  void setIneqPredictionStepStrategy(int strategy);
  void setIneqPredictionParam(double value);

  void setKeepAuxPositiveStrategy(int strategy);
  void setAuxScalingFactor(double value);
  void setAuxCorrectionValue(double value);

  int stepSizeStrategy() const;
  double ineqBacktrackingStepMin() const;
  double auxBacktrackingStepMin() const;

  int initKappaStrategy() const;
  double kappaInitial() const;

  int updateKappaStrategy() const;
  double kappaUpdateFactor() const;
  double kappaFinal() const;
  double tau() const;
  bool limitKappaFinal() const;

  int initSlackStrategy() const;
  double slackVariableInitialIneq() const;

  int initLagrangeStrategy() const;
  double lagrangeMultiplierInitialIneq() const;

  int ineqPredictionStrategy() const;
  int ineqPredictionStepStrategy() const;
  double ineqPredictionParam() const;

  int keepAuxPositiveStrategy() const;
  double auxScalingFactor() const;
  double auxCorrectionValue() const;

protected:
  std::unordered_map<void *, std::function<void(int)>>
      edgeProcessingFunctionMap;

  std::unordered_map<void *, std::function<double()>> getDualityGapFunctionMap;
  // solver parameter
  double
      _kappa; // Default value for the primal duality function parameter _kappa
public:
  // clang-format off

// Step-size selection strategy:
// 0 -> use min(ineq_step, aux_step)
// 1 -> use average(ineq_step, aux_step)
// 2 -> use ineq_step only
// 3 -> use aux_step only
int _step_size_strategy = 0;

// Minimum backtracking step for inequality-related updates.
// Setting this to 1.0 disables backtracking for this part.
double _ineq_backtracking_step_min = 0.5;

// Minimum backtracking step for auxiliary-variable updates.
// Setting this to 1.0 disables backtracking for this part.
double _aux_backtracking_step_min = 0.5;


// Initial barrier parameter strategy:
// 0 -> use constant value _kappa_initial
// 1 -> initialize based on duality gap
int _init_kappa_strategy = 0;

// Initial value of the barrier parameter kappa.
double _kappa_initial = 0.01;


// Barrier update strategy:
// 0 -> increase from previous kappa using tau if duality gap shows insufficient progress
// 1 -> update from initial kappa if duality gap shows insufficient progress
// 2 -> directly accept the value suggested by the duality gap
// 3 -> update using fixed factor _kappa_update_factor (nu in the paper)
int _update_kappa_strategy = 0;

// Fixed update factor for kappa when using the corresponding strategy.
// Typical values are in the range [10, 100].
double _kappa_update_factor = 10.0;

// Target value of the barrier parameter kappa.
double _kappa_final = 1500.0;

// Multiplicative growth factor for kappa when _update_kappa_strategy == 0.
double _tau = 1.2;

// Whether to clamp kappa to _kappa_final or allow higher target.
bool _limit_kappa_final = true;


// Slack initialization strategy:
// 0 -> max(-g(x0), _slack_variable_initial_ineq)
// 1 -> max(abs(g(x0)), _slack_variable_initial_ineq)
// 2 -> constant value _slack_variable_initial_ineq
int _init_slack_strategy = 0;

// Default initial slack value for inequality constraints. Must be > 0.
double _slack_variable_initial_ineq = 10.0;


// Lagrange multiplier initialization strategy:
// 0 -> constant value _lagrange_multiplier_initial_ineq
// 1 -> max(1 / (_kappa_initial * slack), _lagrange_multiplier_initial_ineq)
int _init_lagrange_strategy = 0;

// Default initial Lagrange multiplier value for inequality constraints. Must be > 0.
double _lagrange_multiplier_initial_ineq = 10.0;


// Inequality prediction strategy:
// 0 -> evaluate g(x + zeta_pred * delta_x)
// 1 -> use linearized prediction zeta_pred * J * delta_x + g(x)
int _ineq_prediction_strategy = 0;

// Step-size strategy used in inequality prediction (zeta_pred in the paper):
// 0 -> zeta_pred = 1
// 1 -> zeta_pred = _ineq_prediction_param
// 2 -> zeta_pred = step_size
// 3 -> zeta_pred = step_size * _ineq_prediction_param
int _ineq_prediction_step_strategy = 0;

// Additional scaling parameter used in the inequality prediction rule.
double _ineq_prediction_param = 1.0;


// Strategy to keep auxiliary variables positive:
// 0 -> scale the update
// 1 -> correct values below the threshold to a positive value
int _keep_aux_positive_strategy = 0;

// Scaling factor applied to auxiliary-variable updates when
// positivity preservation uses scaling (beta in the paper).
double _aux_scaling_factor = 0.2;

// Positive threshold used when correcting auxiliary variables (theta in the paper)..
double _aux_correction_value = 1e-3;

  // clang-format on
private:
  int _step_size_index =
      0;                   // index for the step size in the backtracking vector
  double _step_size = 1.0; // step size for the update
  double *_update_step = nullptr; // Update step for the decision variables
};

} // namespace cg2o

#include "sparse_optimizer_ispd.hpp"

#endif // G2O_GRAPH_OPTIMIZER_PDIS_H
