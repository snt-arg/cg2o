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
#include <iostream>
#include <unordered_map>

namespace cg2o {

class TerminationCriteria {
public:
  enum class ConvergenceCriterion {
    UpdateNorm = 0, //||update||
    NewtonDecrement = 1, //||b^T * update||
    GradientNorm = 2 // ||b||
  };

  // Public termination thresholds (accessible in optimizer)
  double epsilon_convergence = 1e-2; // Stopping criterion for gradient norm
  double epsilon_constraint =
      1e-2; // Stopping criterion for individual constraint check

private:
  ConvergenceCriterion _convergenceCriterion = ConvergenceCriterion::NewtonDecrement;

  using ConvergenceFunc = std::function<bool(const Eigen::VectorXd &,
                                             const Eigen::VectorXd &, double)>;
  using ConstraintFeasibilityFunc = std::function<bool(double)>;

  std::unordered_map<ConvergenceCriterion, ConvergenceFunc>
      convergenceStrategies;

public:
  TerminationCriteria() {
    // Initialize solution acceptance strategies
    convergenceStrategies[ConvergenceCriterion::UpdateNorm] =
        [this](const Eigen::VectorXd &, const Eigen::VectorXd &updateVec,
               double epsilon) {
          if (epsilon < 0)
            epsilon = std::abs(epsilon) * epsilon_convergence;
          return updateVec.norm() < epsilon;
        };

    convergenceStrategies[ConvergenceCriterion::NewtonDecrement] =
        [this](const Eigen::VectorXd &bVec, const Eigen::VectorXd &updateVec,
               double epsilon) {
          if (epsilon < 0)
            epsilon = std::abs(epsilon) * epsilon_convergence;
          double decrementSquared = updateVec.transpose() * bVec;
          double decrement = std::abs(decrementSquared);
          return decrement < epsilon;
        };

    convergenceStrategies[ConvergenceCriterion::GradientNorm] =
        [this](const Eigen::VectorXd &bVec, const Eigen::VectorXd &,
               double epsilon) {
          if (epsilon < 0)
            epsilon = std::abs(epsilon) * epsilon_convergence;
          return bVec.norm() < epsilon;
        };
  }

  // Public methods to evaluate termination conditions
  bool verifyConvergence(const Eigen::VectorXd &bVec,
                         const Eigen::VectorXd &updateVec, double epsilon) {
    auto it = convergenceStrategies.find(_convergenceCriterion);
    return it != convergenceStrategies.end() &&
           it->second(bVec, updateVec, epsilon);
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
    the termination criterion is based on the Newton decrement, the gradient
    norm, or the update norm where GradientNorm =0, NewtonDecrement =1 ,
    UpdateNorm = 2

    */
    _convergenceCriterion = static_cast<ConvergenceCriterion>(index);

    if (index == 0) {
      std::cout << "The stopping criterion of the update convergence is "
                   "Gradient Norm."
                << std::endl;
    } else if (index == 1) {
      std::cout << "The stopping criterion of the update convergence is Newton "
                   "Decrement."
                << std::endl;
    } else if (index == 2) {
      std::cout
          << "The stopping criterion of the update convergence is Update Norm."
          << std::endl;
    } else {
      throw std::invalid_argument(
          "Invalid stopping criterion. Valid options are 0:Gradient Norm, "
          "1:Newton Decrement, or 2:Update Norm.");
    }
  };
};

} // namespace cg2o


#endif // CG2O_TERMINATION_CRITERIA_H