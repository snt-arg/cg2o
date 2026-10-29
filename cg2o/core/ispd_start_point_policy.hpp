// Templates of ISPDStartPointPolicy that need the complete solver type.
#ifndef CG2O_ISPD_START_POINT_POLICY_HPP
#define CG2O_ISPD_START_POINT_POLICY_HPP

#include "ispd_start_point_policy.h"
#include "sparse_optimizer_ispd.h"
#include <stdexcept>
#include <type_traits>
#include <utility>

namespace cg2o {

namespace internal {
template <typename T> struct always_false : std::false_type {};
// C++17 detection of the constraint edge kinds (no concepts needed)
template <typename T, typename = void> struct is_ineq_edge : std::false_type {};
template <typename T>
struct is_ineq_edge<
    T, std::void_t<decltype(std::declval<T &>().slackVariable()),
                   decltype(std::declval<T &>().lagrangeMultiplier()),
                   decltype(std::declval<T &>().slackUpdate())>>
    : std::true_type {};
template <typename T, typename = void> struct is_eq_edge : std::false_type {};
template <typename T>
struct is_eq_edge<
    T,
    std::void_t<decltype(std::declval<T &>().getVertexLagrangeMultiplierId()),
                decltype(std::declval<T &>().setLagrangeMultiplier(
                    std::declval<T &>().lagrangeMultiplier()))>>
    : std::true_type {};
} // namespace internal

template <typename T>
void ISPDStartPointPolicy::addShiftSequence(SparseOptimizerISPD &solver,
                                            const std::vector<T *> &stages) {
  _shift_members.check(stages); // also rejects a null entry in a short one
  if (stages.size() < 2) {
    return; // nothing to shift
  }
  for (auto *p : stages) {
    if constexpr (std::is_base_of_v<g2o::OptimizableGraph::Vertex, T>) {
      if (solver.vertex(p->id()) != p ||
          p->dimension() != stages.front()->dimension())
        throw std::invalid_argument("shift vertices must belong to this graph "
                                    "and have equal dimensions");
    } else {
      if (!solver.edges().count(p))
        throw std::invalid_argument("shift edge is not in this graph");
    }
  }
  _shift_members.add(stages);

  // stage k takes the value of stage k + n; the last n stages take the value
  // of the last stage (which keeps its own)
  if constexpr (std::is_base_of_v<g2o::OptimizableGraph::Vertex, T>) {
    for (auto *v : stages)
      solver.registerWarmStartVertex(v);
    // primal vertices: x_k <- x_{k+n}, fixed vertices are left untouched
    _shift_functions.push_back([stages](bool primal, int n) {
      if (!primal) {
        return; // the user owns the primal guess
      }
      for (std::size_t k = 0; k + 1 < stages.size(); ++k) {
        if (!stages[k]->fixed()) {
          stages[k]->setEstimate(
              stages[shiftSource(k, n, stages.size())]->estimate());
        }
      }
    });
  } else if constexpr (internal::is_ineq_edge<T>::value) {
    // inequality edges: (lambda_k, s_k, g_k) <- (lambda_{k+n}, s_{k+n},
    // g_{k+n}); g is the reference of the change detector
    _shift_functions.push_back([this, stages](bool, int n) {
      for (std::size_t k = 0; k + 1 < stages.size(); ++k) {
        auto *from = stages[shiftSource(k, n, stages.size())];
        stages[k]->lagrangeMultiplier() = from->lagrangeMultiplier();
        stages[k]->slackVariable() = from->slackVariable();
        shiftRows(stages[k], from,
                  static_cast<int>(from->lagrangeMultiplier().size()));
      }
    });
  } else if constexpr (internal::is_eq_edge<T>::value) {
    // equality edges: nu_k <- nu_{k+n}
    _shift_functions.push_back([stages](bool, int n) {
      for (std::size_t k = 0; k + 1 < stages.size(); ++k) {
        stages[k]->setLagrangeMultiplier(
            stages[shiftSource(k, n, stages.size())]->lagrangeMultiplier());
      }
    });
  } else {
    static_assert(internal::always_false<T>::value,
                  "addWarmStartShiftSequence: T must be a g2o vertex, a "
                  "BaseFixedSizedEdgeIneq or a BaseFixedSizedEdgeEq");
  }
}

} // namespace cg2o

#endif // CG2O_ISPD_START_POINT_POLICY_HPP
