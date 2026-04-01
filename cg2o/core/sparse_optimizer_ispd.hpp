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
#include <algorithm>
#include <cstddef>
#include <cstdlib>
#include <iostream>
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
  double mu = 1.0 / _t;
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
void SparseOptimizerISPD::constructQuadraticFormEq(
    BaseFixedSizedEdgeEq<D, E, VertexTypes...> &edge) {

  auto error = edge.error();
  auto gamma = edge.lagrangeMultiplier();
  auto rho = edge.rho();
  Eigen::Matrix<double, 2 * D, 1> weightedError =
      Eigen::Matrix<double, 2 * D, 1>::Zero();
  Eigen::DiagonalMatrix<double, 2 * D> omega_matrix =
      Eigen::DiagonalMatrix<double, 2 * D>(
          Eigen::Matrix<double, 2 * D, 1>::Zero());
  for (int i = 0; i < D; ++i) {
    omega_matrix.diagonal()[i] = rho[i];
    weightedError[i] =
        -(omega_matrix.diagonal()[i] * error[i] + 0.5 * gamma[i]);
  }

  static const std::size_t _nr_of_vertices =
      sizeof...(VertexTypes) + 1; // add the lagrange multiplier vertex
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

template <int D, typename E, typename... VertexTypes>
bool SparseOptimizerISPD::addEdgeEqImpl(
    BaseFixedSizedEdgeEq<D, E, VertexTypes...> *e) {
  // create new vertex for the lagrangian
  auto vLagrangian = std::make_shared<VertexLagrangeMultiplier<D>>();
  // auto* nu = new VertexLagrangeMultiplier<D>();
  vLagrangian->setId(e->getVertexLagrangeMultiplierId());
  vLagrangian->setInitialValue(_lagrange_multiplier_initial_eq);

  // Add the Nu vertex to the optimizer
  this->addVertex(vLagrangian.get());

  // assign the vertex to the edge
  e->setVertex(e->vertices().size() - 1, vLagrangian.get());

  bool eresult = OptimizableGraph::addEdge(e);
  if (!eresult) {
    std::cerr << "[Error] adding Eq edge to the optimizer" << std::endl;
    return false;
  }

  vLagrangian->setFixed(false);
  e->initializeInformationMatrix();

  // add teh edge to the set of inequality set
  _edgeEqSet.insert(e);
  _vEqLagrangeMultipliers.push_back(vLagrangian);

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

  auto slack_variable = edge->slackVariable().array().abs();
  auto lagrange_multiplier = edge->lagrangeMultiplier().array().abs();
  return slack_variable.matrix().transpose() * lagrange_multiplier.matrix();
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

    this->initializeSlackVariable<D, E, VertexTypes...>(edge);

    auto &multiplier = edge->lagrangeMultiplier();
    auto &slack = edge->slackVariable();

    switch (_init_lagrange_strategy) {
    case 0: {
      multiplier.setConstant(_lagrange_multiplier_initial_ineq);
      break;
    }
    case 1: {
      multiplier = (1.0 / (_t * slack.array()))
                       .cwiseMax(Eigen::Matrix<double, D, 1>::Constant(
                                     _lagrange_multiplier_initial_ineq)
                                     .array())
                       .matrix();
      break;
    }
    }
    // std::cout << "initial slack: " << slack.transpose() << ", multiplier " <<
    // multiplier.transpose() << std::endl;

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

            auto jacobianWorkspace = this->jacobianWorkspace();
            edge->linearizeOplus(jacobianWorkspace);

            edge->computeJacobianUpdateProduct(update, predicted_error);

            predicted_error.array() *= step_size_local;
            predicted_error.array() +=
                Eigen::Map<const Eigen::VectorXd>(edge->errorData(), D).array();
            break;
          }
          case 1: {               // g (x + step * delta x)
            edge->computeError(); //

            // 1- update the vertices of the edge
            updateEdgeVertices(edge, update, step_size_local);

            edge->computeError(); //

            edge->errorData(); // get the data and continue
            predicted_error =
                Eigen::Map<const Eigen::VectorXd>(edge->errorData(), D);

            updateEdgeVertices(edge, update,
                               -step_size_local); // Reverse the update step

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
        1.0 / (_t * slack.array()) +
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


} // namespace cg2o

#endif // G2O_GRAPH_OPTIMIZER_PDIS_HPP
