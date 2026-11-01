// Primal and dual rows of the KKT system of the constrained optimizers. The
// equality multipliers are the dual variables; every other variable is
// primal.
#ifndef CG2O_KKT_STRUCTURE_H
#define CG2O_KKT_STRUCTURE_H

#include <vector>

namespace cg2o {

// A linear solver that treats primal and dual rows differently (e.g. +sigma
// on the primal and -delta on the dual diagonal) implements this. Before
// every optimize(), the constrained optimizers pass the block indices (the
// Hessian indices) of the equality multipliers, when the solver is used
// through a g2o::BlockSolverX.
class KKTStructureAware {
public:
  virtual ~KKTStructureAware() = default;
  virtual void setDualBlocks(std::vector<int> blocks) = 0;
};

} // namespace cg2o

#endif // CG2O_KKT_STRUCTURE_H
