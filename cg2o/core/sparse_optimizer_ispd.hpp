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

#ifndef G2O_GRAPH_OPTIMIZER_PDIS_HPP
#define G2O_GRAPH_OPTIMIZER_PDIS_HPP

#include "sparse_optimizer_for_constraints.h"
#include "sparse_optimizer_ispd.h"
#include "vertex_lagrange_multiplier.h"
#include "warm_start_math.h"
#include "vertex_state_guard.h"
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <optimizable_graph.h>
#include <vector>

namespace cg2o {

// Implementation of the quadratic form construction
template <int D, typename E, typename... VertexTypes>
void SparseOptimizerISPD::constructQuadraticFormIneq(
    BaseFixedSizedEdgeIneq<D, E, VertexTypes...> &edge) {
  // Inverse of the error (if needed)
  auto error = edge.error();
  auto &lambda = edge.lagrangeMultiplier();
  auto &slack = edge.slackVariable();
  double mu = 1.0 / _kappa;
  Eigen::Matrix<double, D, 1> weightedError;
  Eigen::Matrix<double, D, 1> omega;

  Eigen::Matrix<double, D, 1> slack_inv = slack.cwiseInverse();
  omega = (lambda.array() * slack_inv.array()).matrix();
  weightedError = (-mu * slack_inv.array()).matrix();
  weightedError.array() -= lambda.array();
  weightedError.array() -= omega.array() * error.array();

  auto omega_matrix = omega.asDiagonal(); // Convert into a diagonal matrix
  static const std::size_t _nr_of_vertices = sizeof...(VertexTypes);
  edge.constructQuadraticFormNs(omega_matrix, weightedError,
                                std::make_index_sequence<_nr_of_vertices>());
}


template <int D, typename E, typename... VertexTypes>
bool SparseOptimizerISPD::addEdgeIneqImpl(
    BaseFixedSizedEdgeIneq<D, E, VertexTypes...> *e) {
  bool eresult = OptimizableGraph::addEdge(e);
  if (!eresult) {
    std::cerr << "[Error] adding Ineq edge to the optimizer" << std::endl;
    return false;
  }
  this->initializeSlackVariable<D, E, VertexTypes...>(e);

  Eigen::Matrix<double, D, 1> lagrange_multiplier =
      Eigen::Matrix<double, D, 1>::Constant(_lagrange_multiplier_initial_ineq);

  // set the lagrange multiplier and slack variable
  e->setLagrangeMultiplier(lagrange_multiplier);

  // set the ConstructQuadraticFormImpl function
  auto lambda = [this](BaseFixedSizedEdgeIneq<D, E, VertexTypes...> &edge) {
    this->constructQuadraticFormIneq(edge);
  };
  e->setConstructQuadraticFormImpl(lambda);

  // Store the function in a map instead
  edgeProcessingFunctionMap[e] = [this, e](int controller) {
    this->edgeProcessing<D, E, VertexTypes...>(e, controller);
  };

  getDualityGapFunctionMap[e] = [this, e]() {
    return this->getDualityGap<D, E, VertexTypes...>(e);
  };

  // add teh edge to the set of inequality set
  _edgeIneqSet.insert(e);

  return true;
}

template <int D, typename E, typename... VertexTypes, typename EdgeType>
double SparseOptimizerISPD::getDualityGap(EdgeType *edge) {
  // Ensure EdgeType is  BaseFixedSizedEdgeIneq
  static_assert(
      std::is_base_of<BaseFixedSizedEdgeIneq<D, E, VertexTypes...>,
                      EdgeType>::value,
      "updateMultipliers in PDIS can only accept BaseFixedSizedEdgeIneq");

  // first calculate the lagrange multiplier update and slack variable update
  // update_lagrange_multiplier = lagrange_multiplier / slack_variable *[ g(x) +
  // Jacobian(g(x)) * update_step (x)]  +   1/(_t *slack_variable)
  // update_slack_variable =  - [ s + g(x) + Jacobian(g(x)) * update_step (x)]

  const auto &s = edge->slackVariable();
  const auto &lambda = edge->lagrangeMultiplier();
  for (int i = 0; i < D; ++i) {
    if (!warm_start::validPair(s(i), lambda(i), _kappa))
      return std::numeric_limits<double>::infinity();
  }
  return s.dot(lambda);
}

template <int D, typename E, typename... VertexTypes, typename EdgeType>
void SparseOptimizerISPD::edgeProcessing(EdgeType *edge, int controller) {
  // Ensure EdgeType is  BaseFixedSizedEdgeIneq
  static_assert(
      std::is_base_of<BaseFixedSizedEdgeIneq<D, E, VertexTypes...>,
                      EdgeType>::value,
      "updateMultipliers in PDIS can only accept BaseFixedSizedEdgeIneq");

  switch (controller) {
  case 0: { // reset the lagrange multiplier and slack variable
    this->initializeAuxiliaryVariables<D, E, VertexTypes...>(edge);
    return;
  }

  case 4: { // warm start: store (lambda, slack, g) of the current iterate
    edge->computeError();
    edge->storeAuxiliary(_solve_counter);
    return;
  }

  case 5: { // warm start: restore the stored point (if it is from that solve)
    edge->setAuxiliaryWarm(edge->restoreAuxiliary(_warm_point_stamp));
    return;
  }

  case 6: { // warm start: cold initialisation of edges without warm data
    if (!edge->auxiliaryWarm()) {
      this->initializeAuxiliaryVariables<D, E, VertexTypes...>(edge);
    }
    return;
  }

  case 10: { edge->storeStartingAuxiliary(_solve_counter); return; }
  case 11: { edge->restoreStartingAuxiliary(_solve_counter); return; }
  case 15: { // always checked, independently of optional push / centering
    for (int i = 0; i < D; ++i)
      _aux_valid = _aux_valid && warm_start::validPair(
          edge->slackVariable()(i), edge->lagrangeMultiplier()(i), _kappa);
    return;
  }

  case 7: { // warm start: push s and lambda away from zero
    auto &slack = edge->slackVariable();
    auto &multiplier = edge->lagrangeMultiplier();
    const double ps = std::max(_warm_start_push, _warm_start_slack_push);
    const double pl = std::max(_warm_start_push, _warm_start_multiplier_push);
    edge->computeError();
    for (int i = 0; i < D; ++i) {
      if (!warm_start::validPair(slack(i), multiplier(i), _kappa)) {
        // Infeasible-start repair: do not confuse s>0 with g(x)<0.
        slack(i) = std::max(-edge->errorData()[i],
                           std::max(ps, _slack_variable_initial_ineq));
        multiplier(i) = std::max(pl, _lagrange_multiplier_initial_ineq);
        ++_last_reset_rows;
      }
      if (_warm_start_slack_projection == 1)
        slack(i) = std::max(-edge->errorData()[i], std::max(ps, 1e-12));
      slack(i) = std::max(slack(i), ps);
      multiplier(i) = std::max(multiplier(i), pl);
    }
    return;
  }

  case 8: { // warm start: re-centre the pairs into [lo, hi] / kappa
    if (centringActive()) {
      auto &slack = edge->slackVariable();
      auto &multiplier = edge->lagrangeMultiplier();
      const double mu = 1.0 / _kappa;
      for (int i = 0; i < D; ++i) {
        if (!warm_start::centerPair(slack(i), multiplier(i), mu,
                _warm_start_center_lo, _warm_start_center_hi,
                _warm_start_center_mode,
                std::max(_warm_start_push, _warm_start_slack_push),
                std::max(_warm_start_push, _warm_start_multiplier_push)))
          _aux_valid = false;
      }
    }
    return;
  }

  case 13: { // statistics: products s_i lambda_i of this edge
    const auto &slack = edge->slackVariable();
    const auto &multiplier = edge->lagrangeMultiplier();
    for (int i = 0; i < D; ++i) {
      const double p = slack(i) * multiplier(i);
      _acc_sum += p;
      _acc_min = std::min(_acc_min, p);
      _acc_max = std::max(_acc_max, p);
      ++_acc_rows;
    }
    return;
  }

  case 14: { // strategy 3: the current values are the warm point
    edge->computeError();
    edge->warmIneq() =
        Eigen::Map<const Eigen::Matrix<double, D, 1>>(edge->errorData());
    edge->setAuxiliaryWarm(true);
    return;
  }

  case 9: { // warm start: re-initialise rows whose constraint data changed
    // g of the warm row now (new data, warm primal) is compared with g of
    // the stored iterate (shifted along with lambda and s). Without a data
    // change the two are equal: the primal is shifted with the rows. A large
    // change (e.g. an obstacle moved onto the path, g from -1e4 to > 0)
    // means the stored (lambda, s) of that row describe another constraint:
    // lambda / s ~ 0 would give it almost no weight in the first steps.
    // Such rows are re-initialised at the current x according to
    // _warm_start_reset_mode.
    if (!edge->auxiliaryWarm()) {
      return; // already initialised cold
    }
    edge->computeError();
    auto &slack = edge->slackVariable();
    auto &multiplier = edge->lagrangeMultiplier();
    const auto &gOld = edge->warmIneq();
    const double mu = 1.0 / _kappa;
    const double s_min = std::sqrt(mu);
    for (int i = 0; i < D; ++i) {
      const double g = edge->errorData()[i];
      // test 0: g changed since the stored iterate; test 1: the slack does
      // not match g (primal residual g + s), which also flags rows of an
      // intermediate iterate that were still far from feasible
      const double ref = (_warm_start_reset_test == 0) ? gOld(i) : -slack(i);
      const double scale = _warm_start_reset_test == 2
          ? std::max({std::numeric_limits<double>::min(), std::abs(g), std::abs(slack(i))})
          : std::max({1.0, std::abs(ref), slack(i)});
      if (!(std::abs(g - ref) <= _warm_start_reset_tolerance * scale) ||
          !warm_start::validPair(slack(i), multiplier(i), _kappa)) {
        if (_warm_start_reset_mode == 0) { // centred at mu
          slack(i) = std::max(-g, s_min);
          multiplier(i) = mu / slack(i);
        } else { // cold values of this row
          slack(i) = std::max(-g, _slack_variable_initial_ineq);
          multiplier(i) = _lagrange_multiplier_initial_ineq;
        }
        ++_last_reset_rows;
      }
    }
    return;
  }

  case 1: { // compute the predicted error
    auto computeErrorPrediction =
        [this](EdgeType *edge, Eigen::Matrix<double, D, 1> &predicted_error,
               const double step_size_global, const double *update) {
          double step_size_local;
          switch (_ineq_prediction_step_strategy) {
          case 0: {
            step_size_local = 1.0;
            break;
          }
          case 1: {
            step_size_local = _ineq_prediction_param;
            break;
          }
          case 2: {
            step_size_local = step_size_global;

            break;
          }
          case 3: {
            step_size_local = step_size_global * _ineq_prediction_param;
            break;
          }
          default:
            throw std::runtime_error(
                "[Error] Unknown error prediction step size strategy.");
          }
          switch (_ineq_prediction_strategy) {
          case 0: { // g(x)+ step * Jacobian(g(x)) * delta x
            predicted_error.setZero();
            edge->computeError();

            auto &jacobianWorkspace = this->jacobianWorkspace();
            edge->linearizeOplus(jacobianWorkspace);

            edge->computeJacobianUpdateProduct(update, predicted_error);

            predicted_error.array() *= step_size_local;
            predicted_error.array() +=
                Eigen::Map<const Eigen::VectorXd>(edge->errorData(), D).array();
            break;
          }
          case 1: {               // g (x + step * delta x)
            edge->computeError(); //
            internal::VertexStateGuard trial(edge->vertices());

            // 1- update the vertices of the edge
            updateEdgeVertices(edge, update, step_size_local);

            edge->computeError(); //

            edge->errorData(); // get the data and continue
            predicted_error =
                Eigen::Map<const Eigen::VectorXd>(edge->errorData(), D);

            trial.restore();
            edge->computeError();

            break;
          }
          default:
            throw std::runtime_error(
                "[Error] Unknown error prediction strategy");
          }

          return;
        };

    auto &predicted_error = edge->rho(); // we save the predicted error in rho
    computeErrorPrediction(edge, predicted_error, _step_size,
                           this->_update_step);

    return;
  }
  case 2: { // backtracking for the lagrange multiplier and slack variable

    auto applyBacktracking = [this](const auto &variable, const auto &update) {
      while ((std::size_t)_step_size_index < _alphaBacktracking.size() - 1 &&
             _alphaBacktracking[_step_size_index] >
                 _aux_backtracking_step_min) {
        auto step_size = _alphaBacktracking[_step_size_index];
        if (((variable + step_size * update).array() >= 0).all()) {
          return;
        }
        _step_size_index++;
      }
    };

    auto &slack = edge->slackVariable();
    auto &multipliers = edge->lagrangeMultiplier();
    auto &predicted_error = edge->rho(); // we save the predicted error in rho
    auto &slack_update = edge->slackUpdate();
    auto &multipliers_update = edge->multiplierUpdate();

    slack_update = -slack.array() - predicted_error.array();

    multipliers_update =
        1.0 / (_kappa * slack.array()) +
        multipliers.array() / slack.array() * predicted_error.array();

    applyBacktracking(multipliers, multipliers_update);
    applyBacktracking(slack, slack_update);
    if (false) {
      std::cout << "multipliers : " << multipliers.transpose() << std::endl;
      std::cout << "multipliers update: " << multipliers_update.transpose()
                << std::endl;
      std::cout << "slack : " << slack.transpose() << std::endl;
      std::cout << "slack update: " << slack_update.transpose() << std::endl;
    }
    return;
  };

  case 3: { // update the lagrange multiplier and slack variable
    auto &slack = edge->slackVariable();
    auto &multipliers = edge->lagrangeMultiplier();
    auto &slack_update = edge->slackUpdate();
    auto &multipliers_update = edge->multiplierUpdate();

    // Strategy 1: Scaling factor
    auto updateWithScalingStrategy =
        [this](Eigen::Matrix<double, D, 1> &variable, const auto &update,
               double step_size) {
          for (int i = 0; i < variable.size(); ++i) {
            double new_value = variable(i) + step_size * update(i);
            if (new_value <= 0.0) {
              ++_scaled_count;
              variable(i) = variable(i) * _aux_scaling_factor;
            } else {
              variable(i) = new_value;
            }
          }
        };

    // Strategy 2: Correction to positive value
    auto updateWithCorrectionStrategy =
        [this](Eigen::Matrix<double, D, 1> &variable, const auto &update,
               double step_size) {
          for (int i = 0; i < variable.size(); ++i) {
            double new_value = variable(i) + step_size * update(i);
            if (new_value <= 0.0) {
              ++_scaled_count;
              variable(i) = _aux_correction_value;
              std::cout << "small value" << std::endl;
            } else {
              variable(i) = new_value;
            }
          }
        };

    // Strategy 3: absolute inequality correction
    auto updateWithAbsIneqStrategy =
        [this, edge](Eigen::Matrix<double, D, 1> &variable, const auto &update,
                     double step_size) {
          for (int i = 0; i < variable.size(); ++i) {
            double new_value = variable(i) + step_size * update(i);
            if (new_value <= 0.0) {
              ++_scaled_count;
              variable(i) = std::max(_slack_variable_initial_ineq,
                                     std::abs(edge->errorData()[i]));
              std::cout << "small value" << std::endl;
            } else {
              variable(i) = new_value;
            }
          }
        };
    /****************************************************************** */

    // Choose which strategy to use
    switch (_keep_aux_positive_strategy) {
    case 0:
      updateWithScalingStrategy(multipliers, multipliers_update, _step_size);
      updateWithScalingStrategy(slack, slack_update, _step_size);
      break;
    case 1:
      updateWithCorrectionStrategy(multipliers, multipliers_update, _step_size);
      updateWithCorrectionStrategy(slack, slack_update, _step_size);
      break;
    case 2:
      updateWithAbsIneqStrategy(multipliers, multipliers_update, _step_size);
      updateWithAbsIneqStrategy(slack, slack_update, _step_size);
      break;
    default:
      throw std::runtime_error("[Error] Unknown keep_aux_positive_strategy. 0 "
                               "for scaling, 1 for correction.");
    }

    return;
  }
  default:
    throw std::runtime_error("[Error] Unknown controller in edgeProcessing");
  }
}

template <int D, typename E, typename... VertexTypes>
void SparseOptimizerISPD::initializeSlackVariable(
    BaseFixedSizedEdgeIneq<D, E, VertexTypes...> *e) {
  Eigen::Matrix<double, D, 1> slack_variable;

  switch (_init_slack_strategy) {
  case 0: {
    // set the slack variable to max(-g(x0),_slack_variable_initial_ineq)
    e->computeError();
    for (int i = 0; i < D; ++i) {
      double g_x0 = -(e->errorData()[i]);
      if (g_x0 > _slack_variable_initial_ineq) {
        slack_variable[i] = g_x0;
      } else {
        slack_variable[i] = _slack_variable_initial_ineq;
      }
    }
    break;
  }
  case 1: {
    // set the slack variable to
    // max(abs(g(x0)),_slack_variable_initial_ineq)
    e->computeError();
    for (int i = 0; i < D; ++i) {
      double g_x0 = std::abs(e->errorData()[i]);
      if (g_x0 > _slack_variable_initial_ineq) {
        slack_variable[i] = g_x0;
      } else {
        slack_variable[i] = _slack_variable_initial_ineq;
      }
    }
    break;
  }
  case 2: {
    // set the slack variable to constant _slack_variable_initial_ineq
    slack_variable =
        Eigen::Matrix<double, D, 1>::Constant(_slack_variable_initial_ineq);
    break;
  }
  default:
    throw std::runtime_error("[Error] Unknown slack initialization strategy. "
                             "0: max(-g(x0),_slack_variable_initial_ineq)"
                             "1: max(abs(g(x0)),_slack_variable_initial_ineq)"
                             "2: constant  =  _slack_variable_initial_ineq");
  }

  e->setSlackVariable(slack_variable);
};




template <int D, typename E, typename... VertexTypes>
void SparseOptimizerISPD::initializeAuxiliaryVariables(
    BaseFixedSizedEdgeIneq<D, E, VertexTypes...> *edge) {
  this->initializeSlackVariable<D, E, VertexTypes...>(edge);

  auto &multiplier = edge->lagrangeMultiplier();
  auto &slack = edge->slackVariable();

  switch (_init_lagrange_strategy) {
  case 0: {
    multiplier.setConstant(_lagrange_multiplier_initial_ineq);
    break;
  }
  case 1: {
    multiplier = (1.0 / (_kappa * slack.array()))
                     .cwiseMax(Eigen::Matrix<double, D, 1>::Constant(
                                   _lagrange_multiplier_initial_ineq)
                                   .array())
                     .matrix();
    break;
  }
  }
}

namespace internal {
template <typename T> struct always_false : std::false_type {};
// C++17 detection of the constraint edge kinds (no concepts needed)
template <typename T, typename = void>
struct is_ineq_edge : std::false_type {};
template <typename T>
struct is_ineq_edge<
    T, std::void_t<decltype(std::declval<T &>().slackVariable()),
                   decltype(std::declval<T &>().lagrangeMultiplier()),
                   decltype(std::declval<T &>().warmIneq()),
                   decltype(std::declval<T &>().setAuxiliaryWarm(true))>>
    : std::true_type {};
template <typename T, typename = void> struct is_eq_edge : std::false_type {};
template <typename T>
struct is_eq_edge<
    T, std::void_t<decltype(std::declval<T &>().getVertexLagrangeMultiplierId()),
                   decltype(std::declval<T &>().setLagrangeMultiplier(
                       std::declval<T &>().lagrangeMultiplier()))>>
    : std::true_type {};
} // namespace internal

template <typename T>
void SparseOptimizerISPD::registerWarmStartVertex(T* v) {
  if (!v || vertex(v->id()) != v)
    throw std::invalid_argument("warm-start vertex is not in this graph");
  _snapshot_factories[v] = [v]() {
    using Estimate = std::decay_t<decltype(v->estimate())>;
    Estimate saved = v->estimate();
    return std::function<bool()>([v, saved = std::move(saved)]() {
      v->setEstimate(saved);
      return true;
    });
  };
}

template <typename T>
void SparseOptimizerISPD::addWarmStartShiftSequence(
    const std::vector<T *> &stages) {
  if (stages.size() < 2) {
    return; // nothing to shift
  }
  for (auto *p : stages) {
    if (p == nullptr) {
      throw std::invalid_argument(
          "[Error] addWarmStartShiftSequence: null entry in the sequence");
    }
    if (_shift_members.count(p))
      throw std::invalid_argument("warm-start shift sequences overlap");
  }
  std::unordered_set<T*> unique(stages.begin(), stages.end());
  if (unique.size() != stages.size())
    throw std::invalid_argument("duplicate entry in warm-start shift sequence");
  for (auto* p : stages) {
    if constexpr (std::is_base_of_v<g2o::OptimizableGraph::Vertex, T>) {
      if (vertex(p->id()) != p || p->dimension() != stages.front()->dimension())
        throw std::invalid_argument("shift vertices must belong to this graph and have equal dimensions");
    } else {
      if (!edges().count(p))
        throw std::invalid_argument("shift edge is not in this graph");
    }
  }
  _shift_members.insert(stages.begin(), stages.end());

  // stage k takes the value of stage k + n; the last n stages take the value
  // of the last stage (which keeps its own)
  auto source = [](std::size_t k, std::size_t n, std::size_t size) {
    return std::min(k + n, size - 1);
  };
  if constexpr (std::is_base_of_v<g2o::OptimizableGraph::Vertex, T>) {
    for (auto* v : stages) registerWarmStartVertex(v);
    // primal vertices: x_k <- x_{k+n}, fixed vertices are left untouched
    _shift_functions.push_back([stages, source](bool primal, int n) {
      if (!primal) {
        return; // the user owns the primal guess
      }
      for (std::size_t k = 0; k + 1 < stages.size(); ++k) {
        if (!stages[k]->fixed()) {
          stages[k]->setEstimate(
              stages[source(k, n, stages.size())]->estimate());
        }
      }
    });
  } else if constexpr (internal::is_ineq_edge<T>::value) {
    // inequality edges: (lambda_k, s_k) <- (lambda_{k+n}, s_{k+n})
    _shift_functions.push_back([stages, source](bool, int n) {
      for (std::size_t k = 0; k + 1 < stages.size(); ++k) {
        auto *from = stages[source(k, n, stages.size())];
        stages[k]->lagrangeMultiplier() = from->lagrangeMultiplier();
        stages[k]->slackVariable() = from->slackVariable();
        stages[k]->warmIneq() = from->warmIneq();
        stages[k]->setAuxiliaryWarm(from->auxiliaryWarm());
      }
    });
  } else if constexpr (internal::is_eq_edge<T>::value) {
    // equality edges: nu_k <- nu_{k+n}
    _shift_functions.push_back([stages, source](bool, int n) {
      for (std::size_t k = 0; k + 1 < stages.size(); ++k) {
        stages[k]->setLagrangeMultiplier(
            stages[source(k, n, stages.size())]->lagrangeMultiplier());
      }
    });
  } else {
    static_assert(internal::always_false<T>::value,
                  "addWarmStartShiftSequence: T must be a g2o vertex, a "
                  "BaseFixedSizedEdgeIneq or a BaseFixedSizedEdgeEq");
  }
}

template <int D, typename E, typename... VertexTypes>
bool SparseOptimizerISPD::addEdgeEqImpl(
    BaseFixedSizedEdgeEq<D, E, VertexTypes...> *e) {
  // create new vertex for the lagrangian
  auto *vLagrangian = new VertexLagrangeMultiplier<D>();
  // auto* nu = new VertexLagrangeMultiplier<D>();
  vLagrangian->setId(e->getVertexLagrangeMultiplierId());
  vLagrangian->setInitialValue(_lagrange_multiplier_initial_eq);

  // Add the Nu vertex to the optimizer
  this->addVertex(vLagrangian);

  // assign the vertex to the edge
  e->setVertex(e->vertices().size() - 1, vLagrangian);

  bool eresult = OptimizableGraph::addEdge(e);
  if (!eresult) {
    std::cerr << "[Error] adding Eq edge to the optimizer" << std::endl;
    return false;
  }

  vLagrangian->setFixed(false);
  e->initializeInformationMatrix();
  this->installEqualityQuadraticForm(e); // optional -sigma I block

  // add teh edge to the set of inequality set
  _edgeEqSet.insert(e);
  _vEqLagrangeMultipliers.push_back(vLagrangian);

  return true;
}


} // namespace cg2o

#endif // G2O_GRAPH_OPTIMIZER_PDIS_HPP
