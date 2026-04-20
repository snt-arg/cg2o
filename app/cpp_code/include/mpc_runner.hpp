/*
Copyright (c) 2023, University of Luxembourg
All rights reserved.

Redistributions and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its
   contributors may be used to endorse or promote products derived from
   this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS 'AS IS'
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
*/

#ifndef MPC_RUNNER_H
#define MPC_RUNNER_H

#include "app_config.h"

#include "mpc_formulation.h"
#include "mpc_parameters.h"

// include all inequality optimizers
#include "cg2o/core/sparse_optimizer_al.h"
#include "cg2o/core/sparse_optimizer_bipm.h"
#include "cg2o/core/sparse_optimizer_ispd.h"

#include <deque>
#include <memory>
#include <stdexcept>
#include <type_traits>
#include <variant>

template <class Var, class F> decltype(auto) visit_runner(Var &v, F &&f) {
  return std::visit([&](auto &r) -> decltype(auto) { return f(r); }, v);
}

struct OptimizerType {
  using OptVar = std::variant<std::shared_ptr<g2o::SparseOptimizerAL>,
                              std::shared_ptr<g2o::SparseOptimizerBIPM>,
                              std::shared_ptr<g2o::SparseOptimizerISPD>>;

  OptVar optimizer;

  g2o::SparseOptimizer *base() {
    return std::visit(
        [](auto &opt) -> g2o::SparseOptimizer * { return opt.get(); },
        optimizer);
  }

  template <typename EdgeT> bool addEdgeEq(EdgeT *edge) {
    return std::visit([&](auto &opt) { return opt->addEdgeEq(edge); },
                      optimizer);
  }

  template <typename EdgeT> bool addEdgeIneq(EdgeT *edge) {
    return std::visit([&](auto &opt) { return opt->addEdgeIneq(edge); },
                      optimizer);
  }
};

#endif // MPC_RUNNER_H