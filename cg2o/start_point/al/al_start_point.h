// Start point of SparseOptimizerAL: the multipliers and penalties a solve
// begins with; the primal variables are the values in the graph. The
// optimizer owns one, opt.startPoint().
//
// New AL starts (for example shifted MPC multipliers or least-squares
// multiplier estimates) are added as methods here, without changing the AL
// inner and outer loops. A method that stores and restores an AL point needs
// the multipliers, the per-row penalties and the stored constraint
// violations. The equality multiplier vertices of AL are fixed: they are
// solver state, not measurements, and VertexEstimateCache (which skips fixed
// vertices) does not copy them as it is.
#ifndef CG2O_AL_START_POINT_H
#define CG2O_AL_START_POINT_H

#include "cg2o/start_point/start_point_common.h"

namespace cg2o {

class SparseOptimizerAL;
class ALStartPoint;

//! the built-in AL start methods
enum class ALStartMethod {
  Custom = -1, //!< a method of the application, selected with use<S>()
  Cold = 0,    //!< initial multipliers and penalties
  //! the multipliers and penalties the edges hold, those of the previous
  //! solve (not a stored copy: they are reused in place)
  ReuseMultipliers = 1
};

// What an AL method can do with the start state of the solver.
// ALStartPoint creates one for the duration of a call to the method.
class ALStartContext {
public:
  //! initial multipliers and penalties of the equality and inequality edges
  void resetMultipliers();

private:
  friend class ALStartPoint;
  explicit ALStartContext(ALStartPoint &start) : _start(start) {}
  ALStartPoint &_start;
};

// Interface of the AL start methods.
class ALStartStrategy {
public:
  virtual ~ALStartStrategy() = default;
  //! the built-in method, Custom for others
  virtual ALStartMethod method() const { return ALStartMethod::Custom; }
  //! sets the multipliers and penalties before the iterations
  virtual void prepare(ALStartContext &start) = 0;
};

namespace al {

class Cold final : public ALStartStrategy {
public:
  struct Settings {};
  Settings settings;

  ALStartMethod method() const override { return ALStartMethod::Cold; }
  void prepare(ALStartContext &start) override { start.resetMultipliers(); }
};

class ReuseMultipliers final : public ALStartStrategy {
public:
  struct Settings {};
  Settings settings;

  ALStartMethod method() const override {
    return ALStartMethod::ReuseMultipliers;
  }
  void prepare(ALStartContext &) override {}
};

} // namespace al

class ALStartPoint {
public:
  explicit ALStartPoint(SparseOptimizerAL &solver);
  ALStartPoint(const ALStartPoint &) = delete;
  ALStartPoint &operator=(const ALStartPoint &) = delete;

  //! a built-in method; std::invalid_argument for Custom or an unknown
  //! value, and the method stays
  void select(ALStartMethod method);
  //! makes S the method, also one of the application; returns its settings
  template <typename S> typename S::Settings &use() {
    return _method.template use<S>();
  }
  //! settings of the selected method S; std::logic_error for another one
  template <typename S> typename S::Settings &settings() {
    return _method.template settings<S>();
  }
  ALStartMethod method() const { return _method.active().method(); }

private:
  friend class SparseOptimizerAL;
  friend class ALStartContext;
  //! called by optimize() before the iterations
  void prepare();
  void resetMultipliers();

  SparseOptimizerAL &_solver;
  StrategySlot<ALStartStrategy> _method;
};

} // namespace cg2o

#endif // CG2O_AL_START_POINT_H
