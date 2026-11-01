// Start point of SparseOptimizerBIPM: the barrier parameter and the equality
// multipliers a solve begins with; the primal variables are the values in
// the graph. The optimizer owns one, opt.startPoint().
//
// BIPM needs a strictly feasible primal point; the solver checks it before
// every solve. A method that restores or shifts the primal point must keep
// it strictly feasible for the current problem data: a new kappa cannot
// repair an infeasible point.
#ifndef CG2O_BIPM_START_POINT_H
#define CG2O_BIPM_START_POINT_H

#include "cg2o/start_point/start_point_common.h"

namespace cg2o {

class SparseOptimizerBIPM;
class BIPMStartPoint;

//! the built-in BIPM start methods
enum class BIPMStartMethod {
  Custom = -1, //!< a method of the application, selected with use<S>()
  Cold = 0     //!< kappa_initial, initial equality multipliers
};

//! settings shared by the methods
struct BIPMStartOptions {
  //! a cold start keeps the equality multipliers of the previous solve
  bool keepEqualityMultipliers = false;
};

// What a BIPM method can do with the start state of the solver.
// BIPMStartPoint creates one for the duration of a call to the method.
class BIPMStartContext {
public:
  //! kappa_initial (kappa_final without inequality edges) and the initial
  //! equality multipliers, kept with options().keepEqualityMultipliers
  void coldStart();

private:
  friend class BIPMStartPoint;
  explicit BIPMStartContext(BIPMStartPoint &start) : _start(start) {}
  BIPMStartPoint &_start;
};

// Interface of the BIPM start methods.
class BIPMStartStrategy {
public:
  virtual ~BIPMStartStrategy() = default;
  //! the built-in method, Custom for others
  virtual BIPMStartMethod method() const { return BIPMStartMethod::Custom; }
  //! sets kappa and the equality multipliers before the iterations
  virtual void prepare(BIPMStartContext &start) = 0;
};

namespace bipm {

class Cold final : public BIPMStartStrategy {
public:
  struct Settings {};
  Settings settings;

  BIPMStartMethod method() const override { return BIPMStartMethod::Cold; }
  void prepare(BIPMStartContext &start) override { start.coldStart(); }
};

} // namespace bipm

class BIPMStartPoint {
public:
  explicit BIPMStartPoint(SparseOptimizerBIPM &solver);
  BIPMStartPoint(const BIPMStartPoint &) = delete;
  BIPMStartPoint &operator=(const BIPMStartPoint &) = delete;

  //! a built-in method; std::invalid_argument for Custom or an unknown
  //! value, and the method stays
  void select(BIPMStartMethod method);
  //! makes S the method, also one of the application; returns its settings
  template <typename S> typename S::Settings &use() {
    return _method.template use<S>();
  }
  //! settings of the selected method S; std::logic_error for another one
  template <typename S> typename S::Settings &settings() {
    return _method.template settings<S>();
  }
  BIPMStartMethod method() const { return _method.active().method(); }
  BIPMStartOptions &options() { return _options; }
  const BIPMStartOptions &options() const { return _options; }

private:
  friend class SparseOptimizerBIPM;
  friend class BIPMStartContext;
  //! called by optimize() before the iterations
  void prepare();
  void coldStart();

  SparseOptimizerBIPM &_solver;
  StrategySlot<BIPMStartStrategy> _method;
  BIPMStartOptions _options;
};

} // namespace cg2o

#endif // CG2O_BIPM_START_POINT_H
