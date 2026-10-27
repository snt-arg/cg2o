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
#include <functional>
#include <iostream>
#include <limits>
#include <optimizable_graph.h>
#include <type_traits>
#include <unordered_map>
#include <unordered_set>
#include <utility>
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
   * Karush–Kuhn–Tucker (KKT) conditions.
   *
   * The solver enables joint optimization of:
   * - primal variables (graph states),
   * - Lagrange multipliers,
   * - slack variables,
   * - and barrier parameters,
   * within a unified optimization loop, eliminating the need for nested
   * outer/inner iterations and strictly feasible initialization.
   *
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

  using g2o::SparseOptimizer::addVertex;
  // Typed insertion enables independent snapshots, including manifolds.
  // Base-pointer insertion requires registerWarmStartVertex<T>() or the
  // vertex's estimate-data interface when a warm strategy is selected.
  template <typename T>
    requires std::is_base_of_v<g2o::OptimizableGraph::Vertex, T> &&
             requires(T* v) { v->estimate(); v->setEstimate(v->estimate()); }
  bool addVertex(T* v) {
    if (!g2o::SparseOptimizer::addVertex(v)) return false;
    registerWarmStartVertex(v);
    return true;
  }
  template <typename T> void registerWarmStartVertex(T* v);
  bool removeVertex(g2o::HyperGraph::Vertex* v, bool detach = false) override;
  bool removeEdge(g2o::HyperGraph::Edge* e) override;
  void clear() override;

 
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

  template <int D, typename E, typename... VertexTypes>
  void initializeAuxiliaryVariables(
      BaseFixedSizedEdgeIneq<D, E, VertexTypes...> *edge);

  void executeEdgeProcessing(void *edgePtr, int controller);
  virtual void update(const double *update) override;
  virtual int optimize(int iterations, bool online = false) override;
  void computeKappa(bool input = false);
  double constraintsBacktracking(const double *update);

  // ==========================================================================
  // Warm start
  // ==========================================================================
  //
  // ISPD uses g(x) <= 0, g(x) + s = 0, s > 0, lambda > 0,
  // and s_i*lambda_i approximately 1/kappa. Reuse must respect all of them.
  // 0: cold duals/slacks/barrier, application-provided x.
  // 1: previous final point; usually efficient for small perturbations.
  // 2: previous middle point; a heuristic for avoiding boundary blocking.
  // 3: current values supplied by an application warm-start policy.
  //
  // Owned estimate snapshots coexist with g2o push/pop backups. Insert
  // vertices using their concrete pointer type, or register them explicitly
  // with registerWarmStartVertex<T>(); serialized estimate data is a fallback.
  // Manifold estimates are copied as estimates, never as tangent increments.
  //
  // Register disjoint ordered stage sequences for MPC shifts. Fixed vertices
  // retain the latest measurement. Replicating the last stage is a generic
  // tail heuristic, not a terminal-feasibility guarantee. Use the post-shift
  // callback for a model-consistent terminal rollout. Removing graph objects
  // invalidates the warm cache and clears shift registrations and callbacks.
  //
  // Invalid auxiliary pairs are repaired even if optional pushes, centering,
  // and change detection are off. Failure rollback restores x, nu, s, lambda,
  // and kappa; an optional retry uses cold duals from that restored x. Retry
  // work counts toward optimize()'s returned total (up to twice the budget).
  //
  // No strategy or numerical threshold is universally best. Compare objective,
  // primal feasibility, stationarity, complementarity, tail runtime and retries
  // on the same prescribed problem sequence.
  //
  // Typical receding-horizon use:
  //   opt.setWarmStartStrategy(1);
  //   opt.setWarmStartShift(true);
  //   opt.addWarmStartShiftSequence(states);
  //   opt.addWarmStartShiftSequence(inputs);
  //   opt.addWarmStartShiftSequence(dynamicsEdges);
  //   opt.addWarmStartShiftSequence(inputBoundEdges);
  //   opt.setWarmStartPostShiftCallback(terminalRollout);
  //   // Update fixed x0 before each optimize(). Do not shift x a second time.

  //! 0: cold (default), 1: final iterate, 2: middle iterate,
  //! 3: the values currently in the graph (x, nu, lambda, s and the current
  //!    kappa): for warm-start policies written by the application (load
  //!    any point, then optimize); push, centring, kappa strategy, change
  //!    detector, shift and the breakdown safeguards still apply
  void setWarmStartStrategy(int strategy);
  //! shift the restored point one stage (MPC); needs registered sequences
  void setWarmStartShift(bool shift);
  //! stages per shift (default 1): stage k takes stage k + n, the last n
  //! stages take the last one. Use n > 1 when one control period spans n
  //! stages (e.g. a space-domain horizon at speed: n ~ v dt / ds).
  void setWarmStartShiftSteps(int n);
  int warmStartShiftSteps() const;
  //! restore (and shift) the primal variables too (default true)
  void setWarmStartPrimal(bool primal);
  //! which iterate strategy 2 stores (always the last one that qualifies;
  //! the first iterate is stored if none qualifies):
  //! 0: kappa used in the iteration <= middle kappa
  //! 1: measured duality-gap kappa m / sum(s_i lambda_i) <= middle kappa
  //! 2: iteration index < middle iteration
  //! 3 (default): relative duality gap s'lambda / (1 + |F|) >= middle gap
  //!    (F = robust objective / 2). Inspired by approximate-center LP warm
  //!    starts; 3e-2 is a tunable heuristic. This quantity is affected by
  //!    objective offsets and scaling and is not an NLP optimality certificate.
  void setWarmStartMiddleStrategy(int strategy);
  void setWarmStartMiddleGap(double relativeGap);
  //! optional centrality preference for the middle point: among the
  //! qualifying iterates, the last one inside the neighbourhood N(lo, hi) =
  //! { lo <= s_i lambda_i / mu <= hi for all rows } is stored; if none is
  //! inside, the last qualifying one (Gondzio & Grothey store points in
  //! N(0.5, 2)). lo = 0, hi = inf (default): no preference.
  void setWarmStartMiddleNeighbourhood(double lo, double hi);
  //! middle strategies 0, 1: middle kappa (default 30)
  void setWarmStartMiddleKappa(double kappa);
  void setWarmStartMiddleIteration(int iteration);
  //! kappa of the warm start (clamped to [kappa_initial, kappa_final]):
  //! 0: kappa stored with the point * factor (default)
  //! 1: duality gap of the warm point, factor * m / sum(s_i lambda_i): the
  //!    kappa for which the point is on the central path on average
  //!    (s_i lambda_i = 1/kappa)
  //! 2: constant warm-start kappa value
  //! 3: relative-gap target: kappa = m / (gap * (1 + |F|)), i.e. mu0 =
  //!    gap (1 + |F|) / m at the warm primal. This is a heuristic adaptation
  //!    of LP relative-gap initialization, not a general NLP guarantee.
  void setWarmStartKappaStrategy(int strategy);
  void setWarmStartKappaFactor(double factor);
  void setWarmStartKappaValue(double kappa);
  void setWarmStartKappaGap(double relativeGap);
  //! lower bound for s and lambda of the warm point (0: off, default)
  void setWarmStartPush(double push);
  //! Independent floors, in the units of s and lambda respectively.
  void setWarmStartSlackPush(double push);
  void setWarmStartMultiplierPush(double push);
  //! 0: retain slack; 1: s = max(-g(x), slack push, 1e-12).
  //! Projection reduces the slack residual, but changes complementarity.
  void setWarmStartSlackProjection(int strategy);
  //! centrality band beta >= 1 for s_i lambda_i kappa (0: off, default);
  //! same as setWarmStartCenteringBand(1 / beta, beta)
  void setWarmStartCentering(double beta);
  //! re-centre rows with s_i lambda_i kappa outside [lo, hi] onto the
  //! nearest band edge (lo = hi = 1: exactly centred; lo = 0, hi = inf: off)
  void setWarmStartCenteringBand(double lo, double hi);
  //! how a row is moved into the band:
  //! 0 (default): s_i and lambda_i are scaled by the same factor, which
  //!    keeps s_i / lambda_i unless explicit component floors intervene
  //! 1: legacy s_i < lambda_i rule; depends on units (not an active-set test)
  //! 2: preserve s_i, adjust lambda_i (useful with slack projection)
  //! 3: preserve lambda_i, adjust s_i
  //! Floors take precedence if the requested band and floors conflict.
  void setWarmStartCenteringMode(int mode);
  //! true (default): after a solve that did not converge (iteration cap or
  //! linear-solver failure) the next solve starts cold
  void setWarmStartFallback(bool fallback);
  //! true (default): if a warm-started solve fails (linear-solver failure
  //! after divergence), optimize() puts the starting point back and solves
  //! again with a cold start of nu, lambda, s and kappa (worst case twice the
  //! iterations). With false the starting point is put back and optimize()
  //! returns 0, as for any failed solve.
  void setWarmStartRetry(bool retry);
  //! a solve that does not converge and ends with a constraint violation
  //! above factor * max(violation at its start, epsilon_constraint) counts as
  //! broken down, like a linear-solver failure: the starting point is put
  //! back (and retried cold if it was a warm start). 0: only failures count.
  //! Default 10. Applies whenever a warm-start strategy (1, 2) is selected.
  void setWarmStartDivergenceFactor(double factor);
  double warmStartDivergenceFactor() const;
  //! change detector: a warm row whose g at the start of the solve differs
  //! from g at the stored iterate (both at the same, shifted, primal) by
  //! more than tol * max(1, |g_stored|, s) is treated as a changed
  //! constraint (see setWarmStartResetMode). 0: off. Default 0.5.
  void setWarmStartResetTolerance(double tolerance);
  //! which reference the change detector compares g with:
  //! 0: g of the stored iterate (only a change of the constraint data),
  //! 1 (default): -s, the primal residual |g + s|. For a converged (final)
  //!    point s = -g, so both tests flag the same rows. For a middle point
  //!    test 1 also flags rows that were not yet consistent, and the solve
  //!    then starts with cold duals from the middle primal; in the tests
  //!    this made the middle strategy robust to sudden obstacles.
  //! 2: |g+s| / max(|g|, |s|, smallest normal double): a relative residual
  //!    without the unit offset. Invariant under positive constraint-row
  //!    rescaling, away from the numerical lower bound. Useful for small bounds.
  void setWarmStartResetTest(int test);
  int warmStartResetTest() const;
  //! 0: flagged rows centred at mu, 1: flagged rows get cold values,
  //! 2 (default): any flagged row -> cold start of nu, lambda, s and kappa
  //! (primal stays warm)
  void setWarmStartResetMode(int mode);
  int warmStartResetMode() const;

  int warmStartStrategy() const;
  bool warmStartShift() const;
  bool warmStartPrimal() const;
  int warmStartMiddleStrategy() const;
  double warmStartMiddleKappa() const;
  int warmStartMiddleIteration() const;
  int warmStartKappaStrategy() const;
  double warmStartKappaFactor() const;
  double warmStartKappaValue() const;
  double warmStartPush() const;
  double warmStartCentering() const;
  bool warmStartFallback() const;
  bool warmStartRetry() const;
  double warmStartResetTolerance() const;

  //! registers one stage sequence for the MPC shift (stage 0 first)
  template <typename T>
  void addWarmStartShiftSequence(const std::vector<T *> &stages);
  void clearWarmStartShiftSequences();
  //! After all shifts, before repairs: e.g. roll out the terminal state with
  //! the terminal input. Cleared with the shift sequences on graph removal.
  void setWarmStartPostShiftCallback(std::function<void()> callback) {
    _post_shift_callback = std::move(callback);
  }
  //! forget the stored point; the next solve starts cold
  void resetWarmStart();
  //! true if a stored point exists that the next optimize() can use
  bool warmStartAvailable() const;
  //! what the last optimize() started from: 0 cold, 1 final, 2 middle
  int lastWarmStartUsed() const;
  //! barrier parameter at the end of the last optimize()
  double lastKappa() const;
  //! kappa the last optimize() started with
  double lastStartKappa() const;
  //! true if the last optimize() met the convergence test
  bool lastConverged() const;
  //! true if the last optimize() had to retry cold after a failed warm solve
  bool lastRetried() const;
  //! Actual iterations across both attempts, also when optimize() returns 0
  //! to indicate failure (the historical return-value convention).
  int lastIterationCount() const { return _last_iteration_count; }
  //! rows re-initialised by the reset tolerance in the last optimize()
  int lastResetRows() const;

  //! true if the last optimize() broke down and its starting point was put
  //! back into the vertices (the vertices do not hold the last iterate)
  bool lastRestored() const;

  // ---- complementarity / centrality statistics --------------------------
  // For the starting point and for each iterate: mu = s'lambda / m, the
  // spread of the products s_i lambda_i around mu (centrality C = max_i
  // |s_i lambda_i / mu - 1|, as in the long-step neighbourhood N(lo, hi)),
  // the relative gap s'lambda / (1 + |cost|), and the step sizes and
  // safeguard activity of the iteration.
  struct PointStats {
    int iteration = -1;     //!< -1: starting point of the solve
    int rows = 0;           //!< inequality rows m
    double kappa = 0;       //!< barrier parameter used (target 1/kappa)
    double mu = 0;          //!< s'lambda / m
    double minRatio = 0;    //!< min_i s_i lambda_i / mu
    double maxRatio = 0;    //!< max_i s_i lambda_i / mu
    double centrality = 0;  //!< max_i |s_i lambda_i / mu - 1|
    double relativeGap = std::numeric_limits<double>::quiet_NaN(); //!< s'lambda / (1+|F|)
    double cost = std::numeric_limits<double>::quiet_NaN(); //!< F = chi2/2 of the linearisation point
    double violation = std::numeric_limits<double>::quiet_NaN(); //!< max(|h|, g+), start only
    double ineqStep = 1.0;  //!< inequality backtracking step of the iteration
    double auxStep = 1.0;   //!< lambda/s backtracking step of the iteration
    double appliedStep = 1.0; //!< step applied to x
    int attempt = 0;        //!< 0 initial solve, 1 cold retry
    int scaled = 0;         //!< lambda/s components hit by the positivity safeguard
  };
  //! starting point of the last optimize() (computed when a warm strategy
  //! is selected or the statistics are recorded)
  const PointStats &lastStartStats() const;
  //! first iteration of the last optimize() (first step sizes)
  const PointStats &lastFirstIterationStats() const;
  //! record the statistics of every iterate (off by default; O(m) each)
  void setRecordIterationStats(bool record);
  const std::vector<PointStats> &iterationStats() const;

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
      _kappa = 0.0; // initialized on every cold solve
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
// 0 -> use linearized prediction zeta_pred * J * delta_x + g(x)
// 1 -> evaluate g(x + zeta_pred * delta_x) (nonlinear prediction heuristic)
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


// ---- Warm start (see setWarmStartStrategy) ----
// 0 -> cold start, 1 -> final iterate, 2 -> middle iterate of the previous solve
int _warm_start_strategy = 0;
// Shift registered stage sequences before the solve (MPC), by n stages.
bool _warm_start_shift = false;
int _warm_start_shift_steps = 1;
// Restore (and shift) the primal variables as well.
bool _warm_start_primal = true;
// Middle iterate: 0 -> kappa <= _warm_start_middle_kappa,
// 1 -> measured duality-gap kappa <= _warm_start_middle_kappa,
// 2 -> iteration index < _warm_start_middle_iteration,
// 3 -> relative duality gap >= _warm_start_middle_gap (default)
int _warm_start_middle_strategy = 3;
double _warm_start_middle_kappa = 30.0;
int _warm_start_middle_iteration = 2;
double _warm_start_middle_gap = 3e-2;
double _warm_start_middle_lo = 0.0; // neighbourhood filter (off)
double _warm_start_middle_hi = std::numeric_limits<double>::infinity();
// Warm kappa: 0 -> stored kappa * factor, 1 -> factor * m / gap, 2 -> value
int _warm_start_kappa_strategy = 0;
double _warm_start_kappa_factor = 1.0;
double _warm_start_kappa_value = 100.0;
double _warm_start_kappa_gap = 1e-3;
// Lower bound for s and lambda after restoring (0 = off).
double _warm_start_push = 0.0;
double _warm_start_slack_push = 0.0;
double _warm_start_multiplier_push = 0.0;
int _warm_start_slack_projection = 0;
// Centring band [lo, hi] for s_i lambda_i kappa (lo = 0, hi = inf: off) and
// mode (0 ratio-preserving, 1 keep lambda of active / s of inactive rows).
double _warm_start_center_lo = 0.0;
double _warm_start_center_hi = std::numeric_limits<double>::infinity();
int _warm_start_center_mode = 0;
// Start cold after a solve that did not converge.
bool _warm_start_fallback = true;
// Retry cold (from the restored starting point) when a warm solve fails.
bool _warm_start_retry = true;
// Not converged and violation > factor * max(start violation, eps): broken.
double _warm_start_divergence_factor = 10.0;
// Changed-constraint test: |g_now - ref| > tol * max(1, |ref|, s) with
// ref = g_stored (test 0) or ref = -s (test 1).
double _warm_start_reset_tolerance = 0.5;
int _warm_start_reset_test = 1;
// What happens to such rows: 0 -> centred (s = max(-g, sqrt(mu)), lambda =
// mu / s), 1 -> cold values of the row, 2 -> if any row is flagged, cold
// start of all lambda, s, nu and kappa (the warm primal is kept)
int _warm_start_reset_mode = 2;

  // clang-format on
protected:
  // warm-start helpers (sparse_optimizer_ispd.cpp)
  bool prepareStartingPoint();
  void coldStart(bool forceEqualityReset = false);
  void storeStartingPoint();
  void restoreStartingPoint();
  bool auxiliariesValid();
  void validateWarmStartSettings() const;
  int iterate(int iterations, bool online, bool &warm_stored, bool &failed);
  void storeWarmPoint(double kappa);
  void discardVertexSnapshot(); // both snapshots
  double dualityGapKappa(); // m / sum(s_i lambda_i) over active ineq edges
  PointStats pointStats(bool withCost); // products of the current iterate
  bool centringActive() const;
  bool inMiddleNeighbourhood(const PointStats &st) const;
  // accumulators filled by edge controller 13
  double _acc_sum = 0, _acc_min = 0, _acc_max = 0;
  int _acc_rows = 0;
  bool _aux_valid = true;
  // step sizes and safeguard count of the last update()
  double _last_ineq_step = 1, _last_aux_step = 1, _last_applied_step = 1;
  int _scaled_count = 0;
  PointStats _start_stats, _first_stats;
  bool _record_stats = false;
  std::vector<PointStats> _iteration_stats;
  double constraintViolation(); // max(|h|, g+) over active constraint edges
  bool isEqMultiplierVertex(g2o::OptimizableGraph::Vertex *v);

  struct VertexSnapshotEntry {
    g2o::OptimizableGraph::Vertex *vertex;
    int id;
    bool multiplier; // equality multiplier vertex (nu), not primal
    std::function<bool()> restore;
  };
  void pushSnapshot(std::vector<VertexSnapshotEntry> &list, bool primal);
  void restoreSnapshot(std::vector<VertexSnapshotEntry> &list, bool primal);
  void discardSnapshot(std::vector<VertexSnapshotEntry> &list);
  std::vector<VertexSnapshotEntry> _warm_snapshot;
  std::vector<VertexSnapshotEntry> _start_snapshot;
  using SnapshotFactory = std::function<std::function<bool()>()>;
  std::unordered_map<g2o::OptimizableGraph::Vertex*, SnapshotFactory> _snapshot_factories;
  double _start_point_kappa = 0.0;
  std::vector<std::function<void(bool, int)>> _shift_functions; // primal, n
  std::unordered_set<void*> _shift_members;
  std::function<void()> _post_shift_callback;
  std::unordered_set<g2o::OptimizableGraph::Vertex *> _eq_multiplier_set;
  long _solve_counter = 0;       // id of the current / last solve
  long _warm_point_stamp = -1;   // solve that stored the warm point
  double _warm_point_kappa = 0;  // kappa stored with the point
  int _warm_point_source = 0;    // 1 final, 2 middle
  bool _warm_point_valid = false;
  bool _last_converged = false;
  int _last_warm_used = 0;
  double _last_start_kappa = 0;
  bool _last_retry = false;
  int _last_iteration_count = 0;
  bool _last_restored = false;
  int _last_reset_rows = 0;

private:
  int _step_size_index =
      0;                   // index for the step size in the backtracking vector
  double _step_size = 1.0; // step size for the update
  double *_update_step = nullptr; // Update step for the decision variables
};

} // namespace cg2o

#include "sparse_optimizer_ispd.hpp"

#endif // G2O_GRAPH_OPTIMIZER_PDIS_H
