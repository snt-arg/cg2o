#ifndef MPC_RUNNER_H
#define MPC_RUNNER_H
 
#include "app_config.h"

#include "mpc_formulation.h"
#include "mpc_parameters.h"

// include all inequality optimizers
#include "cg2o/core/sparse_optimizer_al.h"
#include "cg2o/core/sparse_optimizer_bipm.h"
#include "cg2o/core/sparse_optimizer_ispd.h"

#include <variant>
#include <memory>
#include <deque>
#include <type_traits>
#include <stdexcept>

template <class Var, class F>
decltype(auto) visit_runner(Var& v, F&& f) {
  return std::visit([&](auto& r) -> decltype(auto) { return f(r); }, v);
}

struct OptimizerType {
  using OptVar = std::variant<
    std::shared_ptr<g2o::SparseOptimizerAL>,
    std::shared_ptr<g2o::SparseOptimizerBIPM>,
    std::shared_ptr<g2o::SparseOptimizerISPD>
  >;

  OptVar optimizer;

  g2o::SparseOptimizer* base() {
    return std::visit([](auto& opt) -> g2o::SparseOptimizer* {
      return opt.get();
    }, optimizer);
  }

  template <typename EdgeT>
  bool addEdgeEq(EdgeT* edge) {
    return std::visit([&](auto& opt) {
      return opt->addEdgeEq(edge);
    }, optimizer);
  }

  template <typename EdgeT>
  bool addEdgeIneq(EdgeT* edge) {
    return std::visit([&](auto& opt) {
      return opt->addEdgeIneq(edge);
    }, optimizer);
  }
};

#endif // MPC_RUNNER_H