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

#ifndef HELPER_OPTIMIZER_SETUP_H
#define HELPER_OPTIMIZER_SETUP_H

#include <g2o/core/block_solver.h>
#include <g2o/solvers/dense/linear_solver_dense.h>
#include <g2o/solvers/eigen/linear_solver_eigen.h>
#include <g2o/solvers/pcg/linear_solver_pcg.h>

#include <g2o/core/optimization_algorithm_dogleg.h>
#include <g2o/core/optimization_algorithm_gauss_newton.h>
#include <g2o/core/optimization_algorithm_levenberg.h>

#include <iostream>
#include <numeric>
#include <vector>

#ifdef MPC_USE_G2O_SOLVERS
#include <g2o/solvers/cholmod/linear_solver_cholmod.h>
#include <g2o/solvers/csparse/linear_solver_csparse.h>
#endif

#ifdef MPC_USE_SOLVERS_UMFPACK
#include "cg2o/solvers/umfpack/linear_solver_eigen_umfpack_lu.h"
#endif

#ifdef MPC_USE_SOLVERS_PARDISO
#include "cg2o/solvers/pardiso/linear_solver_eigen_pardiso_lu.h"
#endif

#ifdef MPC_USE_SOLVERS_EIGEN
#include "cg2o/solvers/eigen/linear_solver_eigen_bicgstab.h"
#include "cg2o/solvers/eigen/linear_solver_eigen_dense_full_lu.h"
#include "cg2o/solvers/eigen/linear_solver_eigen_dense_partial_lu.h"
#include "cg2o/solvers/eigen/linear_solver_eigen_dense_qr.h"
#include "cg2o/solvers/eigen/linear_solver_eigen_ldlt.h"
#include "cg2o/solvers/eigen/linear_solver_eigen_lscg.h"
#include "cg2o/solvers/eigen/linear_solver_eigen_lu.h"
#include "cg2o/solvers/eigen/linear_solver_eigen_qr.h"
#include "cg2o/solvers/eigen/linear_solver_eigen_svd.h"
#endif

void print_resutls(std::vector<double> &results, int N) {
  Eigen::Map<Eigen::VectorXd> resutls(results.data(), results.size());

  // Print table
  std::cout << "┌─────────┬─────────────────────────────────\n";
  std::cout << "│  Name   │  Value                           \n";
  std::cout << "├─────────┼─────────────────────────────────\n";
  std::cout << "│   v     │   " << resutls.segment(0, N + 1).transpose()
            << "  \n";
  std::cout << "│  d_h    │   " << resutls.segment(N + 1, N + 1).transpose()
            << "  \n";
  std::cout << "│  f_t    │   " << resutls.segment(2 * (N + 1), N).transpose()
            << "  \n";
  std::cout << "│  f_b    │   "
            << resutls.segment(3 * (N + 1) - 1, N).transpose() << "  \n";
  std::cout << "│ slack_1 │   "
            << resutls.segment(4 * (N + 1) - 2, N).transpose() << "  \n";
  std::cout << "│ slack_2 │   "
            << resutls.segment(5 * (N + 1) - 3, N).transpose() << "  \n";
  std::cout << "│ slack_3 │   "
            << resutls.segment(6 * (N + 1) - 4, N).transpose() << "  \n";
  std::cout << "└─────────┴──────────────────────────────────\n";
}
void setupLinearSolver(
    int linearSolverType,
    std::unique_ptr<g2o::LinearSolver<g2o::BlockSolverX::PoseMatrixType>>
        &linearSolver) {
  // Set up the linear solver
  switch (linearSolverType) {
#ifdef MPC_USE_SOLVERS_UMFPACK
  case 0:
    std::cout << "Linear solver:        [0]     (EigenUmfPackLU)" << std::endl;
    linearSolver = std::make_unique<
        cg2o::LinearSolverEigenUmfPackLU<g2o::BlockSolverX::PoseMatrixType>>();
    break;
#endif
#ifdef MPC_USE_SOLVERS_PARDISO
  case 1:
    std::cout << "Linear solver:        [1]     (EigenPardisoLU)" << std::endl;
    linearSolver = std::make_unique<
        cg2o::LinearSolverEigenPardisoLU<g2o::BlockSolverX::PoseMatrixType>>();
    break;
#endif
#ifdef MPC_USE_SOLVERS_EIGEN
  case 2:
    std::cout << "Linear solver:        [2]      (EigenLU)  " << std::endl;
    linearSolver = std::make_unique<
        cg2o::LinearSolverEigenLU<g2o::BlockSolverX::PoseMatrixType>>();
    break;
  case 3:
    std::cout << "Linear solver:        [3]      (EigenQR)  " << std::endl;
    linearSolver = std::make_unique<
        cg2o::LinearSolverEigenLDLT<g2o::BlockSolverX::PoseMatrixType>>();
    break;
  case 4:
    std::cout << "Linear solver:        [4]      (EigenQR)  " << std::endl;
    linearSolver = std::make_unique<
        cg2o::LinearSolverEigenQR<g2o::BlockSolverX::PoseMatrixType>>();
    break;
  case 5:
    std::cout << "Linear solver:        [5]     (EigenDensePartialLU)"
              << std::endl;
    linearSolver = std::make_unique<cg2o::LinearSolverEigenDensePartialLU<
        g2o::BlockSolverX::PoseMatrixType>>();
    break;
  case 6:
    std::cout << "Linear solver:        [6]     (EigenDenseFullLU) "
              << std::endl;
    linearSolver = std::make_unique<cg2o::LinearSolverEigenDenseFullLU<
        g2o::BlockSolverX::PoseMatrixType>>();
    break;
  case 7:
    std::cout << "Linear solver:        [7]     (EigenDenseQR)  " << std::endl;
    linearSolver = std::make_unique<
        cg2o::LinearSolverEigenDenseQR<g2o::BlockSolverX::PoseMatrixType>>();
    break;
  case 8:
    std::cout << "Linear solver:        [8]      (EigenSVD)  " << std::endl;
    linearSolver = std::make_unique<
        cg2o::LinearSolverEigenSVD<g2o::BlockSolverX::PoseMatrixType>>();
    break;
  case 9:
    std::cout << "Linear solver:        [9]      (BiCGSTAB)  " << std::endl;
    linearSolver = std::make_unique<
        cg2o::LinearSolverEigenBiCGSTAB<g2o::BlockSolverX::PoseMatrixType>>();
    break;
  case 10:
    std::cout << "Linear solver:        [10]      (BiCGSTAB)  " << std::endl;
    linearSolver = std::make_unique<
        cg2o::LinearSolverEigenLSCG<g2o::BlockSolverX::PoseMatrixType>>();
    break;
#endif
#ifdef MPC_USE_G2O_SOLVERS
  case 11:
    std::cout << "Linear solver:        [11]      (CSparse)  " << std::endl;
    linearSolver = std::make_unique<
        g2o::LinearSolverCSparse<g2o::BlockSolverX::PoseMatrixType>>();
    break;
  case 12:
    std::cout << "Linear solver:        [12]      (Cholmod)  " << std::endl;
    linearSolver = std::make_unique<
        g2o::LinearSolverCholmod<g2o::BlockSolverX::PoseMatrixType>>();
    break;
#endif
  default:
    throw std::runtime_error("Invalid linear solver type");
    return;
  }
}

void setupAlgorithm(const std::string &solverType,
                    std::unique_ptr<g2o::BlockSolverX> &blockSolver,
                    std::unique_ptr<g2o::OptimizationAlgorithm> &algorithm) {
  // Set up the optimization algorithm
  if (solverType == "gn" || solverType == "GN") {
    algorithm = std::make_unique<g2o::OptimizationAlgorithmGaussNewton>(
        std::move(blockSolver));
  } else if (solverType == "lv" || solverType == "LV") {
    algorithm = std::make_unique<g2o::OptimizationAlgorithmLevenberg>(
        std::move(blockSolver));
  } else if (solverType == "dl" || solverType == "DL") {
    algorithm = std::make_unique<g2o::OptimizationAlgorithmDogleg>(
        std::move(blockSolver));
  } else {
    std::cerr << "Invalid solver type" << std::endl;
    return;
  }
  return;
}

void oneTimeSimulationG2O(auto &param, auto &mpc, auto &optimizer,
                          int numberOfIterations) {
  std::vector<double> vp_prediction = {
      1.85,  2.1625, 2.475, 2.475, 2.475, 2.475, 2.475, 2.475, 2.475, 2.475,
      2.475, 2.475,  2.475, 2.475, 2.475, 2.475, 2.475, 2.475, 2.475, 2.475,
      2.475, 2.475,  2.475, 2.475, 2.475, 2.475, 2.475, 2.475, 2.475, 2.475};

  param->set_vp_prediction(vp_prediction);

  param->set_v_h_0(1.42001);
  param->set_d_h_0(5.8451);
  param->set_f_t_prev(1601.85);
  param->set_f_b_prev(0.0);
  param->set_delta_s(1.42001);
  mpc->setInitialGuess(false); // Set the initial guess
  mpc->setFixedInitialGuess(10.0);
  auto start_time = std::chrono::high_resolution_clock::now();
  // Retrieve and print results
  std::vector<double> results = mpc->getResults();
  std::cout << "Initial conditions:" << std::endl;
  print_resutls(results, param->get_N());

  auto num_iterations = optimizer->optimize(numberOfIterations);
  results = mpc->getResults();
  print_resutls(results, param->get_N());
  auto u_input = mpc->getForceInput(0);
  auto end_time = std::chrono::high_resolution_clock::now();
  auto elapsed_msec = std::chrono::duration_cast<std::chrono::milliseconds>(
      end_time - start_time);

  // cout the number of iterations
  std::cout << "Number of iterations: " << num_iterations << std::endl;
  std::cout << "Optimization done." << std::endl;
  std::cout << "u_input: " << u_input << std::endl;

  std::cout << "Elapsed time: " << elapsed_msec.count() << "ms\n";
}

double g2o_compute_input(auto &mpc, auto &optimizer, int numberOfIterations,
                         double &compute_time_ms, double &num_iterations,
                         double &cost_level) {
#ifdef MPC_FEASIBLE_INITIALIZATION
  // Set a feasible initial guess for the optimization
  mpc->setInitialGuess(false);
#else
  mpc->setFixedInitialGuess(10.0); // Set the initial guess
                                   //    mpc->setInitialGuess(false);
#endif
  //  print_resutls(results, mpcHorizon);
  auto start_time = std::chrono::high_resolution_clock::now();
  num_iterations = optimizer->optimize(numberOfIterations);
  auto end_time = std::chrono::high_resolution_clock::now();

  compute_time_ms =
      std::chrono::duration<double, std::milli>(end_time - start_time).count();
  cost_level = mpc->computeCost();

  double u_input = mpc->getForceInput(0);
  return u_input;
}

double compute_a_p_predicted(auto &a_p_memory) {
  std::vector<double> weight_accel_memory = {0, 0, 0, 0,    0,
                                             0, 0, 0, 0.25, 0.75};
  double predicted_accel =
      std::inner_product(weight_accel_memory.begin(), weight_accel_memory.end(),
                         a_p_memory.begin(), 0.0);
  return predicted_accel;
}

double min_of_vector(const std::vector<double> &data) {
  if (data.empty())
    return std::numeric_limits<double>::quiet_NaN();
  return *std::min_element(data.begin(), data.end());
}

void update_controller_parameters(int k, auto &mpc, auto &param,
                                  auto &a_p_memory, double u_input) {

  // if u_input is nan give -222 and give warrning
  if (std::isnan(u_input)) {
    u_input = -222;
    std::cerr << "[Warning] u_input is nan" << std::endl;
  }

  if (u_input >= 0) {
    param->set_f_t_prev(u_input);
    param->set_f_b_prev(0.0);
  } else {
    param->set_f_t_prev(0.0);
    param->set_f_b_prev(-u_input);
  }

  double v_p = param->get_driving_cycle(k);

  if (k > 0) {
    mpc->simulatorVelocity(u_input); // Update the velocity
    mpc->sensor(v_p, a_p_memory);    // Measure the distance to the preceding
                                     // compute a_p and update the memory
  }
  std::vector<double> vp_prediction = mpc->predictor(
      v_p,
      a_p_memory); // Predict the acceleration of the preceding vehicle
  param->set_vp_prediction(vp_prediction);

  bool use_space_domain =
      param->isSpaceDomain() && (min_of_vector(param->get_vp_prediction()) >
                                 param->get_domain_switch_threshold());

  param->set_use_space_domain(use_space_domain);
  if (use_space_domain) {
    std::cout << "*************************************************************"
                 "*************************************************************"
                 "******** SPACE DOMAIN "
              << std::endl;
  }

  if (param->isSpaceDomain()) {
    v_p = v_p + 0.0;
    double delta_s = 1.0 * param->get_v_h_0() * param->get_delta_t();
    param->set_delta_s(delta_s);
    double delta_s_min = 1.0;
    if (delta_s < delta_s_min) {
      param->set_delta_s(delta_s_min);
    }
  }

  if (k == 0) {
    double v_h = 0.0;
    double d_h = 5.5;
    double f_t_prev = 0.0;
    double f_b_prev = 0.0;
    param->set_v_h_0(v_h);
    param->set_d_h_0(d_h);
    param->set_f_t_prev(f_t_prev);
    param->set_f_b_prev(f_b_prev);
    param->set_v_h_prev(v_h);
    param->set_v_p_prev(v_p);
  }
}

double median(const std::vector<double> &sorted, size_t start, size_t end) {
  size_t len = end - start;
  if (len == 0)
    return std::numeric_limits<double>::quiet_NaN();
  size_t mid = start + len / 2;
  if (len % 2 == 0) {
    return (sorted[mid - 1] + sorted[mid]) / 2.0;
  } else {
    return sorted[mid];
  }
}

/**
 * Returns the mean of the data after removing outliers outside
 * the interval [Q1 - 1.5*IQR, Q3 + 1.5*IQR].
 * If the filtered data is empty, returns NaN.
 */
double meanWithoutOutliers(const std::vector<double> &data) {
  if (data.empty())
    return std::numeric_limits<double>::quiet_NaN();

  // Work on a sorted copy
  std::vector<double> sorted = data;
  std::sort(sorted.begin(), sorted.end());

  size_t n = sorted.size();

  // Compute Q1 and Q3 using the median of lower/upper halves
  size_t half = n / 2; // integer division
  double Q1, Q3;

  if (half == 0) {
    return data[0]; // Only one element, return it as the mean
  }

  if (n % 2 == 0) {
    // Even number of elements: lower half = indices [0, half-1], upper = [half,
    // n-1]
    Q1 = median(sorted, 0, half);
    Q3 = median(sorted, half, n);
  } else {
    // Odd number: lower half = indices [0, half-1], upper = [half+1, n-1]
    Q1 = median(sorted, 0, half);
    Q3 = median(sorted, half + 1, n);
  }

  double IQR = Q3 - Q1;
  double lowerBound = Q1 - 1.5 * IQR;
  double upperBound = Q3 + 1.5 * IQR;

  // Filter the original data (preserve order, but order doesn't matter for
  // mean)
  std::vector<double> filtered;
  filtered.reserve(data.size());
  for (double x : data) {
    if (x >= lowerBound && x <= upperBound) {
      filtered.push_back(x);
    }
  }

  if (filtered.empty()) {
    return std::numeric_limits<double>::quiet_NaN();
  }

  double sum = std::accumulate(filtered.begin(), filtered.end(), 0.0);
  return sum / filtered.size();
}

#endif // HELPER_OPTIMIZER_SETUP_H