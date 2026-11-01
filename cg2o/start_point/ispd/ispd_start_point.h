// Start point of SparseOptimizerISPD: from which numerical state a solve
// begins, what is kept for the next solve, and how a broken-down solve is
// undone. The optimizer owns one, opt.startPoint(), and calls it before and
// after the iterations; graph structure stays with initializeOptimization().
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
// need registerVertex<T>(), or fall back to the serialized estimate data.
//
// MPC shift: register disjoint ordered stage sequences with
// addShiftSequence(). Fixed vertices keep the latest measurement. The last
// stages replicate the last one, a generic tail heuristic and not a
// terminal-feasibility guarantee: use the post-shift callback for a
// model-consistent terminal rollout. Removing graph objects invalidates the
// stored point and clears the shift registrations and the callback.
//
// Typical receding-horizon use:
//   auto &start = opt.startPoint();
//   start.select(ISPDStartMethod::PreviousFinal);
//   start.options().shift.enabled = true;
//   start.addShiftSequence(states);
//   start.addShiftSequence(inputs);
//   start.addShiftSequence(dynamicsEdges);
//   start.addShiftSequence(inputBoundEdges);
//   start.setPostShiftCallback(terminalRollout);
//   // Update fixed x0 before each optimize(). Do not shift x a second time.
// Settings of a method: start.settings<ispd::PreviousMiddle>().rule = ...
//
// No method or threshold is best in general. Compare objective, primal
// feasibility, stationarity, complementarity, tail runtime and retries on
// the same problem sequence.
#ifndef CG2O_ISPD_START_POINT_H
#define CG2O_ISPD_START_POINT_H

#include "cg2o/start_point/ispd/ispd_strategies.h"
#include "cg2o/start_point/start_point_common.h"
#include <algorithm>
#include <functional>
#include <limits>
#include <typeindex>
#include <unordered_map>
#include <vector>

namespace cg2o {

class SparseOptimizerISPD;
struct IneqBlock;

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

// Settings shared by the warm methods, and keepEqualityMultipliers for every
// cold start. Set directly:
//   opt.startPoint().options().shift.enabled = true;
// They are validated when optimize() starts (std::invalid_argument).
struct ISPDStartOptions {
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
    //! a solve (any method but Cold) that does not converge and ends with a
    //! violation above factor * max(violation at its start,
    //! epsilon_constraint) counts as broken down, like a linear-solver
    //! failure. 0: only failures count.
    double divergenceFactor = 10.0;
  } recovery;
};

//! what the last optimize() started from
struct ISPDStartReport {
  //! the method whose point the solve started from; Cold after any cold
  //! start (cold method, nothing stored, fallback, cold duals of the change
  //! detector, a point that could not be repaired)
  ISPDStartMethod used = ISPDStartMethod::Cold;
  double kappa = 0.0; //!< kappa the solve started with
  //! rows re-initialised by the change detector or repaired
  int resetRows = 0;
  //! a broken-down warm solve was retried cold
  bool retried = false;
  //! the solve broke down and its starting point was put back: the vertices
  //! do not hold the last iterate
  bool restored = false;
  //! the starting point (computed for warm methods or recorded statistics)
  ISPDPointStats start;
  //! the first iteration (first step sizes)
  ISPDPointStats firstIteration;
};

// Owns the numerical start state of SparseOptimizerISPD: the selected method,
// the stored point of the next solve (primal and equality-multiplier
// estimates, s, lambda, the stored g and kappa), the rollback point of the
// current solve, the stage sequences of the MPC shift and the report of the
// last solve.
class ISPDStartPoint {
public:
  explicit ISPDStartPoint(SparseOptimizerISPD &solver);
  ~ISPDStartPoint();
  ISPDStartPoint(const ISPDStartPoint &) = delete;
  ISPDStartPoint &operator=(const ISPDStartPoint &) = delete;

  // ---- method ----
  //! a built-in method; std::invalid_argument for Custom or an unknown
  //! value, and the method stays. The active method keeps its settings.
  void select(ISPDStartMethod method);
  //! makes S the method, also one of the application; returns its settings
  template <typename S> typename S::Settings &use() {
    return _method.template use<S>();
  }
  //! settings of the selected method S; std::logic_error for another one
  template <typename S> typename S::Settings &settings() {
    return _method.template settings<S>();
  }
  template <typename S> const typename S::Settings &settings() const {
    return _method.template settings<S>();
  }
  ISPDStartMethod method() const { return _method.active().method(); }
  //! settings shared by the methods: shift, kappa rule, repairs, recovery
  ISPDStartOptions &options() { return _options; }
  const ISPDStartOptions &options() const { return _options; }
  //! std::invalid_argument for invalid options or method settings
  void validate() const;

  // ---- graph ----
  //! registers one stage sequence for the MPC shift (stage 0 first): primal
  //! vertices, inequality edges or equality edges
  template <typename T> void addShiftSequence(const std::vector<T *> &stages);
  void clearShiftSequences();
  //! After all shifts, before repairs: e.g. roll out the terminal state with
  //! the terminal input. Cleared with the shift sequences on graph removal.
  void setPostShiftCallback(std::function<void()> callback) {
    _post_shift_callback = std::move(callback);
  }
  //! typed estimate copies for a vertex added through a base pointer
  template <typename T> void registerVertex(T *v);

  // ---- stored point ----
  //! the next optimize() can start from a point stored by the selected
  //! method
  bool available() const;
  //! forget the stored point; the next solve starts cold
  void invalidate();

  //! what the last optimize() started from
  const ISPDStartReport &lastReport() const { return _report; }

private:
  friend class SparseOptimizerISPD;
  friend class ISPDStartContext;

  // ---- used by SparseOptimizerISPD ----
  //! one optimize(): prepares the starting point, runs the iterations,
  //! keeps the point of the next solve and retries a broken-down warm solve
  //! once, cold. run(attempt) runs the iterations; returns what optimize()
  //! returns.
  template <typename Run> int solve(Run &&run);
  //! statistics wanted from iterate `iteration` of the running attempt
  struct StatsRequest {
    bool stats = false;
    bool cost = false;
  };
  StatsRequest statsRequest(int iteration, bool accepted) const;
  //! the statistics of an iterate, whenever the solver computed them
  void afterIteration(const ISPDPointStats &stats, bool accepted);
  //! optimize() stopped before iterating: cold report, no stored point
  void skipSolve();
  //! s and lambda of a new inequality edge
  void initializeNewRows(IneqBlock &block);
  //! graph mutations: the stored point and the shift registrations are
  //! cleared
  void forgetVertex(const g2o::OptimizableGraph::Vertex *v);
  void forgetEdge(const g2o::OptimizableGraph::Edge *e);
  void clear();

  // s, lambda and g of one inequality edge for the stored and the rollback
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
  enum class FinishAction {
    Accept,    //!< the iterate is the result
    RetryCold, //!< the starting point was put back: solve again, cold
    Reject,    //!< the solve broke down; the stored point is dropped
    Error      //!< the algorithm failed to initialize (optimize() returns -1)
  };

  void prepare();
  bool prepareStartingPoint();
  FinishAction finish(const SolveOutcome &outcome);
  void prepareColdRetry();
  bool brokeDown(bool solverFailed);

  // the solver's cold initialization rules
  void initializeCold(bool forceEqualityReset);
  void initializeColdRows(IneqBlock &block);
  void initializeSlack(IneqBlock &block);

  // operations offered to the methods through ISPDStartContext
  void discardCopies();
  bool storedPointUsable() const;
  double restoreStoredPoint();
  double takeGraphPoint();
  void storeCurrentPoint();
  double currentKappa() const;
  bool pairsValid() const;

  RowCache &rowsOf(const g2o::OptimizableGraph::Edge *edge, int dimension);
  bool isWarm(const g2o::OptimizableGraph::Edge *edge) const;
  void shiftRows(const g2o::OptimizableGraph::Edge *to,
                 const g2o::OptimizableGraph::Edge *from, int dimension);
  void storeRollback();
  void restoreRollback();

  SparseOptimizerISPD &_solver;
  StrategySlot<ISPDStartStrategy> _method;
  ISPDStartOptions _options;
  ISPDStartReport _report;
  VertexEstimateCache _vertices;
  std::unordered_map<const g2o::OptimizableGraph::Edge *, RowCache> _rows;
  std::vector<std::function<void(bool, int)>> _shift_functions; // primal, n
  ShiftMembers _shift_members;
  std::function<void()> _post_shift_callback;

  long _solve_counter = 0;      // id of the current / last solve
  long _warm_point_stamp = -1;  // solve that stored the point
  double _warm_point_kappa = 0; // kappa stored with the point
  // method that stored the point
  std::type_index _warm_point_method = typeid(void);
  bool _warm_point_valid = false;
  double _rollback_kappa = 0;

  // state of the current optimize()
  bool _keep_rollback = false;
  bool _started_warm = false;
  bool _retrying = false;
  bool _stored_this_solve = false; // a point was stored during this solve
  double _start_violation = 0.0;
};

template <typename Run> int ISPDStartPoint::solve(Run &&run) {
  prepare();
  SolveOutcome outcome = run(0);
  int done = outcome.iterations;
  FinishAction action = finish(outcome);
  if (action == FinishAction::RetryCold) {
    prepareColdRetry();
    outcome = run(1);
    if (outcome.iterations >= 0)
      done += outcome.iterations;
    action = finish(outcome);
  }
  switch (action) {
  case FinishAction::Error:
    return -1;
  case FinishAction::Reject:
    return outcome.solverFailed ? 0 : done; // 0 reports a solver failure
  default:
    return done;
  }
}

} // namespace cg2o

#endif // CG2O_ISPD_START_POINT_H
