// Start-point policy of SparseOptimizerAL: which multipliers and penalties a
// solve begins with. Graph structure stays with initializeOptimization().
//
// New AL initialization strategies (for example shifted MPC multipliers,
// penalty reuse or least-squares multiplier estimates) are added as a source
// and options here, without changing the AL inner and outer loops.
#ifndef CG2O_AL_START_POINT_POLICY_H
#define CG2O_AL_START_POINT_POLICY_H

namespace cg2o {

class SparseOptimizerAL;

//! where the multipliers and penalties of a solve come from
enum class ALStartSource {
  Cold = 0,               //!< initial multipliers and penalties
  PreviousMultipliers = 1 //!< multipliers and penalties of the previous solve
};

struct ALStartOptions {
  ALStartSource source = ALStartSource::Cold;
};

class ALStartPointPolicy {
public:
  ALStartOptions &options() { return _options; }
  const ALStartOptions &options() const { return _options; }

  //! sets the multipliers and penalties before the iterations; the primal
  //! variables are the values in the graph
  void prepare(SparseOptimizerAL &solver);

private:
  ALStartOptions _options;
};

} // namespace cg2o

#endif // CG2O_AL_START_POINT_POLICY_H
