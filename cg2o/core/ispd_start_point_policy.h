// Start-point policy of SparseOptimizerISPD: from which numerical state a
// solve begins, what is cached for the next solve, and how a broken-down
// solve is undone. Graph structure stays with initializeOptimization().
//
// ISPD uses g(x) <= 0, g(x) + s = 0, s > 0, lambda > 0 and s_i lambda_i
// approximately 1/kappa; a reused point must respect all of them. Invalid
// pairs (s_i, lambda_i) are repaired even if floors, centring and change
// detection are off. A broken-down solve is undone (x, nu, s, lambda and
// kappa are put back) and, optionally, retried from that x with cold nu,
// lambda, s and kappa; the retry counts toward optimize()'s returned total
// (up to twice the budget).
//
// Estimates are copied as owned values (manifolds as estimates, never as
// tangent increments) and coexist with g2o's push/pop backups. Vertices
// inserted with their concrete type are registered automatically; others
// need SparseOptimizerISPD::registerWarmStartVertex<T>(), or fall back to
// the serialized estimate data.
//
// MPC shift: register disjoint ordered stage sequences with
// SparseOptimizerISPD::addWarmStartShiftSequence(). Fixed vertices keep the
// latest measurement. The last stages replicate the last one, a generic tail
// heuristic and not a terminal-feasibility guarantee: use the post-shift
// callback for a model-consistent terminal rollout. Removing graph objects
// invalidates the cached point and clears the shift registrations and the
// callback.
//
// No source or threshold is best in general. Compare objective, primal
// feasibility, stationarity, complementarity, tail runtime and retries on
// the same problem sequence.
#ifndef CG2O_ISPD_START_POINT_POLICY_H
#define CG2O_ISPD_START_POINT_POLICY_H

#include "start_point_common.h"
#include <algorithm>
#include <functional>
#include <limits>
#include <unordered_map>
#include <vector>

namespace cg2o {

class SparseOptimizerISPD;

// For the starting point and for each iterate: mu = s'lambda / m, the
// spread of the products s_i lambda_i around mu (centrality C = max_i
// |s_i lambda_i / mu - 1|, as in the long-step neighbourhood N(lo, hi)),
// the relative gap s'lambda / (1 + |cost|), and the step sizes and
// safeguard activity of the iteration.
struct ISPDPointStats {
  int iteration = -1;    //!< -1: starting point of the solve
  int rows = 0;          //!< inequality rows m
  double kappa = 0;      //!< barrier parameter used (target 1/kappa)
  double mu = 0;         //!< s'lambda / m
  double minRatio = 0;   //!< min_i s_i lambda_i / mu
  double maxRatio = 0;   //!< max_i s_i lambda_i / mu
  double centrality = 0; //!< max_i |s_i lambda_i / mu - 1|
  double relativeGap =
      std::numeric_limits<double>::quiet_NaN(); //!< s'lambda / (1+|F|)
  double cost =
      std::numeric_limits<double>::quiet_NaN(); //!< F = chi2/2 of the
                                                //!< linearisation point
  double violation =
      std::numeric_limits<double>::quiet_NaN(); //!< max(|h|, g+), start only
  double ineqStep = 1.0;    //!< inequality backtracking step of the iteration
  double auxStep = 1.0;     //!< lambda/s backtracking step of the iteration
  double appliedStep = 1.0; //!< step applied to x
  int attempt = 0;          //!< 0 initial solve, 1 cold retry
  int scaled = 0; //!< lambda/s components hit by the positivity safeguard
};

//! where the starting point of a solve comes from. A point stored under
//! PreviousFinal is not used by PreviousMiddle and vice versa: the first
//! solve after a switch starts cold.
enum class ISPDStartSource {
  Cold = 0, //!< cold nu, lambda, s, kappa; application-provided x
  //! final iterate of the previous solve; usually efficient for small
  //! perturbations
  PreviousFinal = 1,
  //! an earlier iterate of the previous solve (see MiddlePointRule); a
  //! heuristic for not starting on the boundary
  PreviousMiddle = 2,
  //! x, nu, lambda, s and the current kappa as they are in the graph, for
  //! warm-start policies written by the application. Floors, centring, the
  //! kappa rule, the change detector, the shift and the breakdown safeguards
  //! still apply.
  CurrentGraph = 3
};

//! which iterate PreviousMiddle stores: the last one that qualifies, the
//! first iterate if none does
enum class MiddlePointRule {
  BarrierKappa = 0,    //!< kappa used in the iteration <= middle.kappa
  DualityGapKappa = 1, //!< measured m / sum(s_i lambda_i) <= middle.kappa
  Iteration = 2,       //!< iteration index < middle.iteration
  //! s'lambda / (1 + |F|) >= middle.gap, F = robust objective / 2. Inspired
  //! by approximate-centre LP warm starts; it depends on objective offsets
  //! and scaling and is not an NLP optimality certificate.
  RelativeGap = 3
};

//! barrier parameter of a warm start, clamped to [kappa_initial,
//! kappa_final]
enum class WarmKappaRule {
  PreviousTimesFactor = 0, //!< kappa stored with the point * kappa.factor
  //! kappa.factor * m / sum(s_i lambda_i): the kappa for which the point is
  //! on the central path on average (s_i lambda_i = 1/kappa)
  DualityGap = 1,
  Constant = 2, //!< kappa.value
  //! m / (kappa.gap * (1 + |F|)), i.e. mu0 = gap (1 + |F|) / m at the warm
  //! primal: a heuristic adaptation of LP relative-gap initialization
  RelativeGap = 3
};

//! how a row is moved into the centring band; floors take precedence if
//! they conflict with the band
enum class CenteringMode {
  //! s and lambda scaled by the same factor, which keeps s / lambda
  PreserveRatio = 0,
  Legacy = 1,            //!< s < lambda rule; depends on units
  PreserveSlack = 2,     //!< adjust lambda (useful with slack projection)
  PreserveMultiplier = 3 //!< adjust s
};

//! reference the change detector compares g at the start of the solve with
enum class ResetTest {
  //! g of the stored iterate: only a change of the constraint data, scaled
  //! by max(1, |ref|, s)
  StoredConstraint = 0,
  //! -s (the residual |g + s|), scaled by max(1, |ref|, s). For a final
  //! point s = -g, so it flags the same rows as StoredConstraint; for a
  //! middle point it also flags rows that were not yet consistent, which
  //! then start with cold duals (in the tests this made PreviousMiddle
  //! robust to sudden obstacles).
  SlackResidual = 1,
  //! |g + s| / max(|g|, |s|): no unit offset, invariant under positive row
  //! scaling; useful for small bounds
  RelativeSlackResidual = 2
};

//! what happens to the rows the change detector flags
enum class ResetMode {
  CentredRows = 0, //!< s = max(-g, sqrt(mu)), lambda = mu / s
  ColdRows = 1,    //!< cold values of the row
  //! any flagged row: cold nu, lambda, s and kappa (x stays warm)
  ColdDuals = 2
};

// All settings of the start-point policy, set directly:
//   opt.startPointPolicy().options().source = ISPDStartSource::PreviousFinal;
// They are validated when optimize() starts (std::invalid_argument),
// whatever the source.
struct ISPDStartOptions {
  ISPDStartSource source = ISPDStartSource::Cold;
  //! restore (and shift) the primal variables too
  bool warmPrimal = true;
  //! a cold start keeps the equality multipliers of the previous solve
  bool keepEqualityMultipliers = false;

  struct Shift {
    bool enabled = false; //!< shift the registered stage sequences (MPC)
    //! stage k takes stage k + steps, the last steps stages the last one;
    //! >= 0. Use steps > 1 when one control period spans several stages
    //! (e.g. a space-domain horizon at speed: steps ~ v dt / ds).
    int steps = 1;
  } shift;

  struct Middle {
    MiddlePointRule rule = MiddlePointRule::RelativeGap;
    double kappa = 30.0; //!< BarrierKappa, DualityGapKappa; > 0
    int iteration = 2;   //!< Iteration; >= 0
    double gap = 3e-2;   //!< RelativeGap, a tunable heuristic; > 0
    //! among the qualifying iterates, prefer the last one inside N(lo, hi) =
    //! { lo <= s_i lambda_i / mu <= hi for all rows } (Gondzio & Grothey
    //! store points in N(0.5, 2)). 0 <= lo <= 1 <= hi; lo = 0, hi = inf:
    //! no preference.
    double neighbourhoodLo = 0.0;
    double neighbourhoodHi = std::numeric_limits<double>::infinity();
  } middle;

  struct Kappa {
    WarmKappaRule rule = WarmKappaRule::PreviousTimesFactor;
    double factor = 1.0;  //!< PreviousTimesFactor, DualityGap; > 0
    double value = 100.0; //!< Constant; > 0
    double gap = 1e-3;    //!< RelativeGap; > 0
  } kappa;

  //! floors of the warm point, in the units of s and lambda (0: off)
  struct Floors {
    double both = 0.0;       //!< lower bound for s and lambda
    double slack = 0.0;      //!< lower bound for s
    double multiplier = 0.0; //!< lower bound for lambda
    //! s = max(-g(x), slack floor, 1e-12): reduces the slack residual, but
    //! changes complementarity
    bool projectSlack = false;
    double slackFloor() const { return std::max(both, slack); }
    double multiplierFloor() const { return std::max(both, multiplier); }
  } push;

  struct Centering {
    //! rows with s_i lambda_i kappa outside [lo, hi] are moved onto the
    //! nearest band edge. 0 <= lo <= 1 <= hi; lo = hi = 1: exactly centred;
    //! lo = 0, hi = inf: off. A symmetric band beta >= 1 is [1/beta, beta].
    double lo = 0.0;
    double hi = std::numeric_limits<double>::infinity();
    CenteringMode mode = CenteringMode::PreserveRatio;
    bool active() const {
      return lo > 0.0 || hi < std::numeric_limits<double>::infinity();
    }
  } centering;

  struct Reset {
    //! change detector: a warm row whose g at the start of the solve differs
    //! from the reference (see ResetTest) by more than tolerance * scale is a
    //! changed constraint, handled by mode. 0: off.
    double tolerance = 0.5;
    ResetTest test = ResetTest::SlackResidual;
    ResetMode mode = ResetMode::ColdDuals;
  } reset;

  struct Recovery {
    //! after a solve that did not converge (iteration cap or linear-solver
    //! failure), the next solve starts cold
    bool fallback = true;
    //! a broken-down warm solve: its starting point is put back and solved
    //! again with cold nu, lambda, s and kappa (worst case twice the
    //! iterations). false: it is put back and optimize() returns 0.
    bool retry = true;
    //! a solve (any source but Cold) that does not converge and ends with a
    //! violation above factor * max(violation at its start,
    //! epsilon_constraint) counts as broken down, like a linear-solver
    //! failure. 0: only failures count.
    double divergenceFactor = 10.0;
  } recovery;
};

// Owns the numerical start state of SparseOptimizerISPD: the cached point of
// the next solve (primal and equality-multiplier estimates, s, lambda, the
// stored g and kappa), the rollback point of the current solve, the stage
// sequences of the MPC shift, and the report of the last solve.
//
// Lifecycle, driven by optimize():
//   prepare() -> [observeIteration() per iterate] -> finish()
//   finish() == RetryCold -> prepareColdRetry() -> iterations -> finish()
class ISPDStartPointPolicy {
public:
  ISPDStartOptions &options() { return _options; }
  const ISPDStartOptions &options() const { return _options; }
  //! throws std::invalid_argument for invalid warm-start settings
  void validate() const;

  // ---- lifecycle ----
  StartResult prepare(SparseOptimizerISPD &solver);
  //! the iterates are observed to select a middle point
  bool observesIterations() const {
    return _options.source == ISPDStartSource::PreviousMiddle;
  }
  //! the middle-point rule needs the objective of each iterate
  bool needsIterationCost() const {
    return _options.middle.rule == MiddlePointRule::RelativeGap;
  }
  void observeIteration(SparseOptimizerISPD &solver,
                        const ISPDPointStats &stats, int iteration);
  FinishAction finish(SparseOptimizerISPD &solver, const SolveOutcome &outcome);
  //! after RetryCold: cold nu, lambda, s, kappa from the restored point
  void prepareColdRetry(SparseOptimizerISPD &solver);
  //! optimize() stopped before iterating: cold report, no cached point
  void skipSolve();
  //! forget the cached point; the next solve starts cold
  void invalidate();

  // ---- cache and registrations ----
  //! a cached point exists that the next optimize() can use: it was stored
  //! for the selected source (PreviousFinal or PreviousMiddle); a point of
  //! the other source is not used, so a switch costs one cold start
  bool available(bool lastConverged) const {
    const bool cachedSource =
        _options.source == ISPDStartSource::PreviousFinal ||
        _options.source == ISPDStartSource::PreviousMiddle;
    return cachedSource && _warm_point_valid &&
           _warm_point_source == static_cast<int>(_options.source) &&
           (lastConverged || !_options.recovery.fallback);
  }
  template <typename T> void registerVertex(T *v) {
    _vertices.registerVertex(v);
  }
  template <typename T>
  void addShiftSequence(SparseOptimizerISPD &solver,
                        const std::vector<T *> &stages);
  void clearShiftSequences();
  void setPostShiftCallback(std::function<void()> callback) {
    _post_shift_callback = std::move(callback);
  }
  //! graph mutations: the cache and the shift registrations are cleared
  void forgetVertex(const g2o::OptimizableGraph::Vertex *v);
  void forgetEdge(const g2o::OptimizableGraph::Edge *e);
  void clear();

  // ---- report of the last optimize() ----
  //! what it started from: 0 cold, 1 final, 2 middle, 3 current values
  int lastSource() const { return _last_source; }
  bool lastRetried() const { return _last_retry; }
  bool lastRestored() const { return _last_restored; }
  int lastResetRows() const { return _last_reset_rows; }

private:
  // s, lambda and g of one inequality edge for the cached and the rollback
  // point; allocated once per edge
  struct RowCache {
    explicit RowCache(int rows) : dimension(rows), values(6 * rows, 0.0) {}
    int dimension;
    std::vector<double> values;
    long cachedStamp = -1;   // solve that stored the cached rows
    long rollbackStamp = -1; // solve that stored the rollback rows
    bool warm = false;       // the edge holds warm values for the coming solve
    double *cachedSlack() { return values.data(); }
    double *cachedMultiplier() { return values.data() + dimension; }
    double *storedConstraint() { return values.data() + 2 * dimension; }
    //! g of the restored point, shifted with s and lambda
    double *referenceConstraint() { return values.data() + 3 * dimension; }
    double *rollbackSlack() { return values.data() + 4 * dimension; }
    double *rollbackMultiplier() { return values.data() + 5 * dimension; }
  };
  RowCache &rowsOf(const g2o::OptimizableGraph::Edge *edge, int dimension);
  bool isWarm(const g2o::OptimizableGraph::Edge *edge) const;
  void shiftRows(const g2o::OptimizableGraph::Edge *to,
                 const g2o::OptimizableGraph::Edge *from, int dimension);

  bool prepareStartingPoint(SparseOptimizerISPD &solver);
  void storeCachedPoint(SparseOptimizerISPD &solver, double kappa);
  void storeRollback(SparseOptimizerISPD &solver);
  void restoreRollback(SparseOptimizerISPD &solver);
  bool brokeDown(SparseOptimizerISPD &solver, bool solverFailed);
  bool inMiddleNeighbourhood(const ISPDPointStats &stats) const;

  ISPDStartOptions _options;
  VertexEstimateCache _vertices;
  std::unordered_map<const g2o::OptimizableGraph::Edge *, RowCache> _rows;
  std::vector<std::function<void(bool, int)>> _shift_functions; // primal, n
  ShiftMembers _shift_members;
  std::function<void()> _post_shift_callback;

  long _solve_counter = 0;      // id of the current / last solve
  long _warm_point_stamp = -1;  // solve that stored the cached point
  double _warm_point_kappa = 0; // kappa stored with the point
  int _warm_point_source = 0;   // 1 final, 2 middle
  bool _warm_point_valid = false;
  double _rollback_kappa = 0;

  // state of the current optimize()
  bool _keep_rollback = false;
  bool _started_warm = false;
  bool _retrying = false;
  bool _stored_this_solve = false; // a point was cached during this solve
  bool _middle_centred = false;    // a centred middle point was cached
  double _start_violation = 0.0;

  int _last_source = 0;
  bool _last_retry = false;
  bool _last_restored = false;
  int _last_reset_rows = 0;
};

} // namespace cg2o

#endif // CG2O_ISPD_START_POINT_POLICY_H
