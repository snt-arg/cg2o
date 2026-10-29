// Start-point policy of SparseOptimizerISPD: from which numerical state a
// solve begins, what is cached for the next solve, and how a broken-down
// solve is undone. Graph structure stays with initializeOptimization().
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

//! where the starting point of a solve comes from
enum class ISPDStartSource {
  Cold = 0,           //!< cold nu, lambda, s, kappa; application-provided x
  PreviousFinal = 1,  //!< final iterate of the previous solve
  PreviousMiddle = 2, //!< a selected earlier iterate of the previous solve
  CurrentGraph = 3    //!< the values currently in the graph (application)
};

//! which iterate PreviousMiddle stores (the last one that qualifies; the
//! first iterate if none does)
enum class MiddlePointRule {
  BarrierKappa = 0,    //!< kappa used in the iteration <= middle kappa
  DualityGapKappa = 1, //!< m / sum(s_i lambda_i) <= middle kappa
  Iteration = 2,       //!< iteration index < middle iteration
  RelativeGap = 3      //!< s'lambda / (1 + |F|) >= middle gap
};

//! barrier parameter of a warm start, clamped to [kappa_initial,
//! kappa_final]
enum class WarmKappaRule {
  PreviousTimesFactor = 0, //!< kappa stored with the point * factor
  DualityGap = 1,          //!< factor * m / sum(s_i lambda_i)
  Constant = 2,            //!< a fixed value
  RelativeGap = 3          //!< m / (gap * (1 + |F|))
};

//! how a pair is moved into the centring band
enum class CenteringMode {
  PreserveRatio = 0,     //!< s and lambda scaled by the same factor
  Legacy = 1,            //!< s < lambda rule; depends on units
  PreserveSlack = 2,     //!< adjust lambda
  PreserveMultiplier = 3 //!< adjust s
};

//! reference the change detector compares g with
enum class ResetTest {
  StoredConstraint = 0,     //!< g of the stored iterate
  SlackResidual = 1,        //!< -s, scaled by max(1, |ref|, s)
  RelativeSlackResidual = 2 //!< -s, scaled by max(|g|, |s|)
};

//! what happens to the rows the change detector flags
enum class ResetMode {
  CentredRows = 0, //!< s = max(-g, sqrt(mu)), lambda = mu / s
  ColdRows = 1,    //!< cold values of the row
  ColdDuals = 2    //!< any flagged row: cold nu, lambda, s, kappa
};

// All settings of the start-point policy. The setWarmStart*() setters of
// SparseOptimizerISPD write these fields and validate each value; fields set
// directly are validated when optimize() starts.
struct ISPDStartOptions {
  ISPDStartSource source = ISPDStartSource::Cold;
  //! restore (and shift) the primal variables too
  bool warmPrimal = true;
  //! a cold start keeps the equality multipliers of the previous solve
  bool keepEqualityMultipliers = false;

  struct Shift {
    bool enabled = false; //!< shift the registered stage sequences (MPC)
    int steps = 1;        //!< stage k takes stage k + steps
  } shift;

  struct Middle {
    MiddlePointRule rule = MiddlePointRule::RelativeGap;
    double kappa = 30.0;
    int iteration = 2;
    double gap = 3e-2;
    //! preference for iterates inside N(lo, hi); lo = 0, hi = inf: none
    double neighbourhoodLo = 0.0;
    double neighbourhoodHi = std::numeric_limits<double>::infinity();
  } middle;

  struct Kappa {
    WarmKappaRule rule = WarmKappaRule::PreviousTimesFactor;
    double factor = 1.0;
    double value = 100.0;
    double gap = 1e-3;
  } kappa;

  struct Floors {
    double both = 0.0;         //!< lower bound for s and lambda (0: off)
    double slack = 0.0;        //!< lower bound for s, in the units of s
    double multiplier = 0.0;   //!< lower bound for lambda
    bool projectSlack = false; //!< s = max(-g(x), slack floor, 1e-12)
    double slackFloor() const { return std::max(both, slack); }
    double multiplierFloor() const { return std::max(both, multiplier); }
  } push;

  struct Centering {
    //! band for s_i lambda_i kappa; lo = 0, hi = inf: off
    double lo = 0.0;
    double hi = std::numeric_limits<double>::infinity();
    CenteringMode mode = CenteringMode::PreserveRatio;
    bool active() const {
      return lo > 0.0 || hi < std::numeric_limits<double>::infinity();
    }
  } centering;

  struct Reset {
    double tolerance = 0.5; //!< 0: change detector off
    ResetTest test = ResetTest::SlackResidual;
    ResetMode mode = ResetMode::ColdDuals;
  } reset;

  struct Recovery {
    bool fallback = true; //!< after a solve that did not converge, start cold
    bool retry = true;    //!< a broken-down warm solve is retried cold
    //! not converged and violation > factor * max(start violation, eps):
    //! broken down. 0: only failures count.
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
  //! a cached point exists that the next optimize() can use
  bool available(bool lastConverged) const {
    return _warm_point_valid && (lastConverged || !_options.recovery.fallback);
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
