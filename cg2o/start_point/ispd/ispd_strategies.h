// Start methods of SparseOptimizerISPD: the strategy interface and the
// built-in methods. A method decides where a solve starts from and what is
// kept for the next solve. ISPDStartPoint does the rest around it: the MPC
// shift, the change detector, the repairs, the kappa rule and the centring
// of a warm point, the rollback of a broken-down solve and its cold retry.
//
// A new method derives from ISPDStartStrategy, has a Settings type and a
// settings member, and is selected with opt.startPoint().use<Method>().
#ifndef CG2O_ISPD_STRATEGIES_H
#define CG2O_ISPD_STRATEGIES_H

#include <limits>
#include <optional>

namespace cg2o {

class ISPDStartPoint;
struct ISPDStartOptions;

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

//! the built-in start methods; the values are the ones of configuration
//! files. A point stored by one method is not used by another: the first
//! solve after a switch starts cold.
enum class ISPDStartMethod {
  Custom = -1, //!< a method of the application, selected with use<S>()
  Cold = 0,    //!< cold nu, lambda, s, kappa; application-provided x
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

//! what a method reads from every accepted iterate
enum class IterateObservation { None, Complementarity, ComplementarityAndCost };

// What a method can do with the start state of the solver. ISPDStartPoint
// creates one for the duration of a call to the method.
class ISPDStartContext {
public:
  //! cold nu (kept with options().keepEqualityMultipliers), lambda, s and
  //! kappa from the solver's initialization rules
  void coldStart();
  //! a point stored by the selected method exists, and the last solve
  //! converged (or options().recovery.fallback is off)
  bool storedPointUsable() const;
  //! puts the stored point back (x only with options().warmPrimal) and
  //! returns its kappa
  double restoreStoredPoint();
  //! the values in the graph are the warm point; returns the current kappa,
  //! kappa_final if none is set
  double takeGraphPoint();
  //! stores the current iterate and kappa as the point of the next solve
  void storeCurrentPoint();
  //! a point was stored during this solve
  bool storedInThisSolve() const;
  //! kappa of the current iterate
  double kappa() const;
  //! every (s, lambda) pair of the current iterate can enter the KKT system
  bool pairsValid() const;
  const ISPDStartOptions &options() const;

private:
  friend class ISPDStartPoint;
  explicit ISPDStartContext(ISPDStartPoint &start) : _start(start) {}
  ISPDStartPoint &_start;
};

// Interface of the ISPD start methods. One object per optimizer; its
// functions are called once per solve or once per iteration, never per row.
class ISPDStartStrategy {
public:
  virtual ~ISPDStartStrategy() = default;
  //! the built-in method, Custom for others
  virtual ISPDStartMethod method() const { return ISPDStartMethod::Custom; }
  //! std::invalid_argument for invalid settings; checked when optimize()
  //! starts
  virtual void validate() const {}
  //! false for a cold method: its solve is not saved for rollback or
  //! retried, and no starting-point statistics are computed for it
  virtual bool warm() const { return true; }
  //! the next optimize() can start from a point this method stored
  virtual bool available(const ISPDStartContext &) const { return false; }
  //! sets the x, nu, s and lambda the solve starts from and returns their
  //! kappa; the shared preparation is then applied. std::nullopt: the
  //! method made a cold start, which is used as it is.
  virtual std::optional<double> load(ISPDStartContext &start) = 0;
  //! statistics the method reads from every accepted iterate
  virtual IterateObservation observation() const {
    return IterateObservation::None;
  }
  virtual void observe(ISPDStartContext &, const ISPDPointStats &) {}
  //! a solve, or its cold retry, begins
  virtual void beginAttempt() {}
  //! the solve is accepted: keep what the next solve starts from
  virtual void retain(ISPDStartContext &) {}
};

//! which iterate PreviousMiddle stores: the last one that qualifies, the
//! first iterate if none does
enum class MiddlePointRule {
  BarrierKappa = 0,    //!< kappa used in the iteration <= kappa
  DualityGapKappa = 1, //!< measured m / sum(s_i lambda_i) <= kappa
  Iteration = 2,       //!< iteration index < iteration
  //! s'lambda / (1 + |F|) >= gap, F = robust objective / 2. Inspired by
  //! approximate-centre LP warm starts; it depends on objective offsets and
  //! scaling and is not an NLP optimality certificate.
  RelativeGap = 3
};

namespace ispd {

// Cold start: the solver's initialization rules, every solve.
class Cold final : public ISPDStartStrategy {
public:
  struct Settings {};
  Settings settings;

  ISPDStartMethod method() const override { return ISPDStartMethod::Cold; }
  bool warm() const override { return false; }
  std::optional<double> load(ISPDStartContext &start) override;
};

// The final iterate of the previous solve.
class PreviousFinal final : public ISPDStartStrategy {
public:
  struct Settings {};
  Settings settings;

  ISPDStartMethod method() const override {
    return ISPDStartMethod::PreviousFinal;
  }
  bool available(const ISPDStartContext &start) const override;
  std::optional<double> load(ISPDStartContext &start) override;
  void retain(ISPDStartContext &start) override;
};

// An earlier iterate of the previous solve, chosen by settings.rule.
class PreviousMiddle final : public ISPDStartStrategy {
public:
  struct Settings {
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
  };
  Settings settings;

  ISPDStartMethod method() const override {
    return ISPDStartMethod::PreviousMiddle;
  }
  void validate() const override;
  bool available(const ISPDStartContext &start) const override;
  std::optional<double> load(ISPDStartContext &start) override;
  IterateObservation observation() const override;
  void observe(ISPDStartContext &start, const ISPDPointStats &stats) override;
  void beginAttempt() override { _centred = false; }
  void retain(ISPDStartContext &start) override;

private:
  bool inNeighbourhood(const ISPDPointStats &stats) const;
  bool _centred = false; // a centred iterate was stored in this attempt
};

// The values currently in the graph.
class CurrentGraph final : public ISPDStartStrategy {
public:
  struct Settings {};
  Settings settings;

  ISPDStartMethod method() const override {
    return ISPDStartMethod::CurrentGraph;
  }
  std::optional<double> load(ISPDStartContext &start) override;
};

} // namespace ispd
} // namespace cg2o

#endif // CG2O_ISPD_STRATEGIES_H
