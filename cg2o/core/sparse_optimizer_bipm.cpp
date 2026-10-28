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

#include "sparse_optimizer_bipm.h"
#include "g2o/core/optimization_algorithm_with_hessian.h"
#include "g2o/core/sparse_optimizer.h"
#include "g2o/stuff/logger.h"
#include "g2o/stuff/timeutil.h"
#include "optimization_algorithm.h"
#include "solver.h"
#include "vertex_state_guard.h"
#include <Eigen/Core>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>

namespace cg2o {
using namespace std;

// Default constructor
SparseOptimizerBIPM::SparseOptimizerBIPM() { _num_inner_iterations_max = 100; }

// Default destructor
SparseOptimizerBIPM::~SparseOptimizerBIPM() = default;

void SparseOptimizerBIPM::resetLagrangeMultiplierEq() {
  for (auto &vertex : _vEqLagrangeMultipliers) {
    vertex->setToOrigin();
  }
}

// Setter solver paramters
void SparseOptimizerBIPM::setKappa(double kappa) { _kappa = kappa; }
void SparseOptimizerBIPM::setKappaInitial(double kappaInitial) {
  _kappa_initial = kappaInitial;
}
void SparseOptimizerBIPM::setKappaUpdateFactor(double kappaUpdateFactor) {
  _kappa_update_factor = kappaUpdateFactor;
}
void SparseOptimizerBIPM::setKappaFinal(double kappaFinal) {
  _kappa_final = kappaFinal;
}

// getter for the solver parameters
double SparseOptimizerBIPM::kappa() const { return _kappa; }
double SparseOptimizerBIPM::kappaInitial() const { return _kappa_initial; }
double SparseOptimizerBIPM::kappaFinal() const { return _kappa_final; }
double SparseOptimizerBIPM::kappaUpdateFactor() const {
  return _kappa_update_factor;
}

double SparseOptimizerBIPM::backtrackingAlgorithm(const double *update) {
  // Recheck all rows for each candidate: nonlinear feasible sets need not
  // remain feasible when a step accepted for an earlier edge is shortened.
  if (_activeEdgesIneq.empty())
    return 1.0;
  // true if all rows of the edge are strictly negative after the trial step;
  // the vertices are restored exactly when the guard goes out of scope
  auto feasibleAt = [&](OptimizableGraph::Edge *edge, double step) {
    internal::VertexStateGuard trial(edge->vertices());
    updateEdgeVertices(edge, update, step);
    edge->computeError();
    for (int i = 0; i < edge->dimension(); ++i)
      if (!(edge->errorData()[i] < 0.0) || !std::isfinite(edge->errorData()[i]))
        return false;
    return true;
  };
  // The edge that rejected the previous step usually rejects the next one
  // too, so it is tested first: a rejected step then costs about one test.
  // The accepted step is unchanged: the first one feasible for all edges.
  const std::size_t numEdges = _activeEdgesIneq.size();
  std::size_t lastFailed = 0;
  double accepted = 0.0;
  for (double step : _alphaBacktracking) {
    if (!feasibleAt(_activeEdgesIneq[lastFailed], step))
      continue;
    bool feasible = true;
    for (std::size_t k = 0; k < numEdges; ++k) {
      if (k == lastFailed)
        continue;
      if (!feasibleAt(_activeEdgesIneq[k], step)) {
        lastFailed = k;
        feasible = false;
        break;
      }
    }
    if (feasible) {
      accepted = step;
      break;
    }
  }
  // Refresh caches once after the search, with the original vertices restored.
  for (auto *edge : _activeEdgesIneq)
    edge->computeError();
  return accepted;
}

bool SparseOptimizerBIPM::strictlyInterior() {
  for (auto *edge : _activeEdgesIneq) {
    edge->computeError();
    for (int i = 0; i < edge->dimension(); ++i)
      if (!(edge->errorData()[i] < 0.0) || !std::isfinite(edge->errorData()[i]))
        return false;
  }
  return true;
}

void SparseOptimizerBIPM::update(const double *update) {

  const double step = backtrackingAlgorithm(update);
  if (step == 0.0) {
    G2O_WARN("BIPM backtracking found no step that keeps the inequality "
             "constraints strictly feasible");
    _barrier_step_failed = true;
    return;
  }
  applyScaledUpdate(update, step);
}

int SparseOptimizerBIPM::optimize(int iterations, bool online) {
  using namespace g2o;
  BatchStatisticsScope statisticsScope;
  if (iterations < 0)
    throw std::invalid_argument("iteration budget must be nonnegative");
  if (iterations == 0 || terminate())
    return 0;
  if (!fixedConstraintsFeasible()) {
    G2O_WARN("constraints on fixed vertices are violated, nothing to "
             "optimize");
    return 0;
  }
  if (!(_kappa_initial > 0) || !std::isfinite(_kappa_initial) ||
      !(_kappa_final >= _kappa_initial) || !std::isfinite(_kappa_final) ||
      !(_kappa_update_factor > 1) || !std::isfinite(_kappa_update_factor))
    throw std::invalid_argument("invalid BIPM barrier settings");
  if (!strictlyInterior()) {
    G2O_WARN("BIPM needs a strictly feasible starting point (all inequality "
             "constraints < 0)");
    return 0;
  }
  _barrier_step_failed = false;

  if (_ivMap.size() == 0) {
    G2O_WARN("0 vertices to optimize, maybe forgot to call "
             "initializeOptimization()");
    return -1;
  }

  int cjIterations = 0;
  double cumTime = 0;
  bool innerLoopStop = false;
  bool outerLoopStop = false;
  bool ok = true;
  _kappa = _kappa_initial; // reset _kappa to the initial value
  if (!_warm_start_lagrange_multiplier_eq_flag) {
    resetLagrangeMultiplierEq();
  }

  ok = _algorithm->init(online);
  if (!ok) {
    G2O_ERROR("Error while initializing");
    return -1;
  }

  _batchStatistics.clear();
  if (_computeBatchStatistics)
    _batchStatistics.resize(iterations);

  OptimizationAlgorithm::SolverResult result = OptimizationAlgorithm::OK;

  if (_activeEdgesIneq.empty()) {
    // set mu to 0 if there are no inequality edges
    this->setKappa(_kappa_final);
  }

  while (!outerLoopStop) {

    int i = 0;
    innerLoopStop = false;
    // BIPM - Inner loop
    while (!innerLoopStop) {
      i++;
      preIteration(cjIterations);
      if (_computeBatchStatistics) {
        G2OBatchStatistics &cstat = _batchStatistics[cjIterations];
        G2OBatchStatistics::setGlobalStats(&cstat);
        cstat.iteration = cjIterations;
        cstat.numEdges = _activeEdges.size();
        cstat.numVertices = _activeVertices.size();
      }

      double ts = get_monotonic_time();
      result = _algorithm->solve(cjIterations, online);
      ok = (result == OptimizationAlgorithm::OK) && !_barrier_step_failed;
      if (_barrier_step_failed)
        result = OptimizationAlgorithm::Fail;

      bool errorComputed = false;
      if (_computeBatchStatistics) {
        computeActiveErrors();
        errorComputed = true;
        _batchStatistics[cjIterations].chi2 = activeRobustChi2();
        _batchStatistics[cjIterations].timeIteration =
            get_monotonic_time() - ts;
      }

      if (verbose()) {
        double dts = get_monotonic_time() - ts;
        cumTime += dts;
        if (!errorComputed)
          computeActiveErrors();
        cerr << "iteration= " << cjIterations
             << "\t chi2= " << FIXED(activeRobustChi2()) << "\t time= " << dts
             << "\t cumTime= " << cumTime
             << "\t edges= " << _activeEdges.size();
        _algorithm->printVerbose(cerr);
        cerr << endl;
      }
      ++cjIterations;
      postIteration(cjIterations);

      // termination criteria
      innerLoopStop = i >= _num_inner_iterations_max ||
                      cjIterations >= iterations || terminate() ||
                      verifyConvergence(-1 * 10) || !ok;
    }

    if (result == OptimizationAlgorithm::Fail) {
      return 0;
    }

    if (cjIterations >= iterations ||
        terminate()) { // check if the maximum number of iterations is reached
      break;
    }
    // check if the termination condition is satisfied
    outerLoopStop =
        verifyConvergence(-1) && verifyEqFeasibility(_activeEdgesEq, -1.0) &&
        verifyIneqFeasibility(_activeEdgesIneq, -1.0); // defualt * 1
    outerLoopStop = _kappa >= _kappa_final && (outerLoopStop);

    // update _kappa
    _kappa = (_kappa < _kappa_final) ? (_kappa * _kappa_update_factor) : _kappa;

  } // End of inner loop

  return cjIterations;
}

} // namespace cg2o
