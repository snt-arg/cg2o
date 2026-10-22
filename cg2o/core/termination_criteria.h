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
#ifndef CG2O_TERMINATION_CRITERIA_H
#define CG2O_TERMINATION_CRITERIA_H


#include "optimizable_graph.h"
#include <Eigen/Dense>
#include <functional>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <unordered_map>

namespace cg2o {

class TerminationCriteria {
public:
  enum class ConvergenceCriterion {
    UpdateNorm = 0, //||update||
    NewtonDecrement = 1, //||b^T * update||
    GradientNorm = 2, // ||b||
    // |b^T * update| < eps * (relative_cost_offset + |cost|): the Newton
    // decrement test scaled by the objective value, so that the tolerance does
    // not depend on the scale of the cost. Needs the current cost; use the
    // four-argument verifyConvergence (the optimizer passes activeChi2()).
    RelativeNewtonDecrement = 3
  };

  // Public termination thresholds (accessible in optimizer)
  double epsilon_convergence = 1e-2; // Stopping criterion for gradient norm
  double epsilon_constraint =
      1e-2; // Stopping criterion for individual constraint check
  // Offset added to |cost| in the RelativeNewtonDecrement test. With a cost
  // near zero the test falls back to the absolute Newton-decrement test.
  double relative_cost_offset = 1.0;

private:
  ConvergenceCriterion _convergenceCriterion = ConvergenceCriterion::NewtonDecrement;

  // arguments: b, update, epsilon, cost
  using ConvergenceFunc = std::function<bool(
      const Eigen::VectorXd &, const Eigen::VectorXd &, double, double)>;
  using ConstraintFeasibilityFunc = std::function<bool(double)>;

  std::unordered_map<ConvergenceCriterion, ConvergenceFunc>
      convergenceStrategies;

public:
  TerminationCriteria() {
    // Initialize solution acceptance strategies
    convergenceStrategies[ConvergenceCriterion::UpdateNorm] =
        [this](const Eigen::VectorXd &, const Eigen::VectorXd &updateVec,
               double epsilon, double) {
          if (epsilon < 0)
            epsilon = std::abs(epsilon) * epsilon_convergence;
          return updateVec.norm() < epsilon;
        };

    convergenceStrategies[ConvergenceCriterion::NewtonDecrement] =
        [this](const Eigen::VectorXd &bVec, const Eigen::VectorXd &updateVec,
               double epsilon, double) {
          if (epsilon < 0)
            epsilon = std::abs(epsilon) * epsilon_convergence;
          double decrementSquared = updateVec.transpose() * bVec;
          double decrement = std::abs(decrementSquared);
          return decrement < epsilon;
        };

    convergenceStrategies[ConvergenceCriterion::GradientNorm] =
        [this](const Eigen::VectorXd &bVec, const Eigen::VectorXd &,
               double epsilon, double) {
          if (epsilon < 0)
            epsilon = std::abs(epsilon) * epsilon_convergence;
          return bVec.norm() < epsilon;
        };

    convergenceStrategies[ConvergenceCriterion::RelativeNewtonDecrement] =
        [this](const Eigen::VectorXd &bVec, const Eigen::VectorXd &updateVec,
               double epsilon, double cost) {
          if (epsilon < 0)
            epsilon = std::abs(epsilon) * epsilon_convergence;
          const double decrement = std::abs(updateVec.dot(bVec));
          return decrement < epsilon * (relative_cost_offset + std::abs(cost));
        };
  }

  // Public methods to evaluate termination conditions
  // cost: current objective value, only used by RelativeNewtonDecrement
  // (cost = 0 gives eps * relative_cost_offset, i.e. the absolute test).
  bool verifyConvergence(const Eigen::VectorXd &bVec,
                         const Eigen::VectorXd &updateVec, double epsilon,
                         double cost) {
    auto it = convergenceStrategies.find(_convergenceCriterion);
    return it != convergenceStrategies.end() &&
           it->second(bVec, updateVec, epsilon, cost);
  }

  // Original signature, unchanged behaviour for the existing criteria.
  bool verifyConvergence(const Eigen::VectorXd &bVec,
                         const Eigen::VectorXd &updateVec, double epsilon) {
    return verifyConvergence(bVec, updateVec, epsilon, 0.0);
  }

  ConvergenceCriterion convergenceCriterion() const {
    return _convergenceCriterion;
  }

  // true when the selected criterion needs the objective value
  bool needsCost() const {
    return _convergenceCriterion ==
           ConvergenceCriterion::RelativeNewtonDecrement;
  }

  bool verifyEqFeasibility(g2o::OptimizableGraph::EdgeContainer &edgeVec,
                           double epsilon) {
    if (epsilon < 0)
      epsilon = std::abs(epsilon) * epsilon_constraint;

    for (auto &edge : edgeVec) {
      edge->computeError();
      const double *error = edge->errorData();
      int edgeDimension = edge->dimension() /
                          2; // divide by 2 to account for equality constraints
                             // excluding the lagrange multiplier
      for (int i = 0; i < edgeDimension; ++i) {
        if (std::abs(error[i]) > epsilon) {
          return false; // If any error exceeds epsilon, return false
        }
      }
    }

    return true; // Return true only if all edges satisfy the condition
  }

  bool verifyIneqFeasibility(g2o::OptimizableGraph::EdgeContainer &edgeVec,
                             double epsilon) {
    if (epsilon < 0)
      epsilon = std::abs(epsilon) * epsilon_constraint;

    for (auto &edge : edgeVec) {
      edge->computeError();
      const double *error = edge->errorData();
      int edgeDimension = edge->dimension();
      for (int i = 0; i < edgeDimension; ++i) {
        if (std::max(0.0, error[i]) > epsilon) {
          return false; // If any error exceeds epsilon, return false
        }
      }
    }

    return true; // Return true only if all edges satisfy the condition
  }

  // Setters to configure termination criteria
  void setConvergenceCriterion(ConvergenceCriterion type) {
    _convergenceCriterion = type;
  }

  void setConvergenceCriterion(int index) {
    /*
    0: UpdateNorm, 1: NewtonDecrement, 2: GradientNorm,
    3: RelativeNewtonDecrement
    */
    const char *name = nullptr;
    switch (index) {
    case 0: name = "Update Norm"; break;
    case 1: name = "Newton Decrement"; break;
    case 2: name = "Gradient Norm"; break;
    case 3: name = "Relative Newton Decrement"; break;
    default:
      throw std::invalid_argument(
          "Invalid stopping criterion. Valid options are 0:Update Norm, "
          "1:Newton Decrement, 2:Gradient Norm, or 3:Relative Newton "
          "Decrement.");
    }
    _convergenceCriterion = static_cast<ConvergenceCriterion>(index);
    std::cout << "The stopping criterion of the update convergence is " << name
              << "." << std::endl;
  }
};

} // namespace cg2o


#endif // CG2O_TERMINATION_CRITERIA_H