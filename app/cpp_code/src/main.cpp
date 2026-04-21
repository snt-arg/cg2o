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

/*
 * Test the ACC MPC
 */

#include "mpc_formulation.h"
#include "mpc_parameters.h"
#include <Eigen/src/Core/GenericPacketMath.h>
#include <Eigen/src/Core/Map.h>
#include <Eigen/src/Core/Matrix.h>
#include <Eigen/src/Core/products/Parallelizer.h>
#include <Eigen/src/Core/util/XprHelper.h>
#include <cmath>
#include <cstdlib>
#include <g2o/core/base_fixed_sized_edge.h>

#include <Eigen/Core>
#include <iostream>
#include <memory>
#include <regex>
#include <string>
#include <unistd.h>
#include <vector>

#include <data_saver.hpp>

#ifdef USE_AL
#include "cg2o/core/sparse_optimizer_al.h"
#endif
#ifdef USE_BIPM
#include "cg2o/core/sparse_optimizer_bipm.h"
#endif
#ifdef USE_ISPD
#include "cg2o/core/sparse_optimizer_ispd.h"
#endif

#include "helper_config.hpp"
#include "helper_optimizer_setup.hpp"

int main(int argc, char **argv) {
  std::cout << "MPC-based ACC" << std::endl;
  std::string config_path = getConfigPath(argc, argv);
  std::cout << "Loading config: " << config_path << std::endl;
  YAML::Node config = YAML::LoadFile(config_path);

  int argCount = 0;

  argCount++;
  int ampl_solver =
      (argc > argCount)
          ? std::atoi(argv[argCount])
          : getIntValue(config["acc_settings"],
                        "ampl_solver"); // 0: ipopt, 1:Couenne, 2: bonmin, 3:
                                        // cplex, 4: gurobi
  // 5: mosek, 6: knitrok, 7: gecode, 8: xpress;

  argCount++;
  int linearSolverType = (argc > argCount)
                             ? std::atoi(argv[argCount])
                             : getIntValue(config["linear_solver"], "type");

  argCount++;
  int mpcHorizon = (argc > argCount) ? std::atoi(argv[argCount])
                                     : getIntValue(config["mpc"], "horizon");

  std::filesystem::path current_path = std::filesystem::current_path();
  // std::cout << "Current path: " << current_path << std::endl;
  // Parse command-line arguments

  std::string solverType =
      getStringValue(config["unconstrained_solver"], "type");
  int numberOfIterations = getIntValue(config["mpc"], "number_of_iterations");
  int run_times = getIntValue(config["mpc"], "run_times");

  std::unique_ptr<g2o::BlockSolverX> blockSolver;
  std::unique_ptr<g2o::LinearSolver<g2o::BlockSolverX::PoseMatrixType>>
      linearSolver;
  setupLinearSolver(linearSolverType, linearSolver);

  blockSolver = std::make_unique<g2o::BlockSolverX>(std::move(linearSolver));
  std::unique_ptr<g2o::OptimizationAlgorithm> algorithm;
  setupAlgorithm(solverType, blockSolver, algorithm);

  // Initialize MPC parameters
  auto param = std::make_shared<cg2o::mpc::MPCParameters>();
  param->set_N(mpcHorizon);
  // using the config to change
  // the parameteres

#ifdef USE_AL
  std::string AlgEqType = "al"; // Change to string, default "al"
  std::string AlIneqType = "AL";
  auto optimizer = std::make_shared<cg2o::SparseOptimizerAL>();
  auto mpc = std::main ampl coloum do you know how data is saved
      ke_shared<cg2o::mpc::MPCFormulation<cg2o::SparseOptimizerAL>>(
          mpcHorizon, optimizer, param);

#endif
#ifdef USE_BIPM
  std::string AlgEqType = "gn"; // Change to string, default "gn"
  std::string AlIneqType = "BIPM";
  auto optimizer = std::make_shared<cg2o::SparseOptimizerBIPM>();
  auto mpc =
      std::make_shared<cg2o::mpc::MPCFormulation<cg2o::SparseOptimizerBIPM>>(
          mpcHorizon, optimizer, param);
#endif
#ifdef USE_ISPD
  std::string AlgEqType = "gn"; // Change to string, default "gn"
  std::string AlIneqType = "ISPD";
  auto optimizer = std::make_shared<cg2o::SparseOptimizerISPD>();
  auto mpc =
      std::make_shared<cg2o::mpc::MPCFormulation<cg2o::SparseOptimizerISPD>>(
          mpcHorizon, optimizer, param);
#endif
  optimizer->setAlgorithm(algorithm.get());
  configureMPCParameters(config, param, optimizer);
  mpc->setupMPC();
  // Update parameters
  // Set the initial guess
  // mpc->setInitialGuess(true);

  std::vector<double> initial = mpc->getResults();
  // std::cout << "Initial Guess:\n"
  //   << Eigen::Map<Eigen::VectorXd>(initial.data(),
  //   initial.size()).transpose() << std::endl;

  int terminationCriterion =
      1; // 0 UpdateNorm   1  NewtonDecrement     2 GradientNorm
  optimizer->setConvergenceCriterion(terminationCriterion);

  optimizer->initializeOptimization();
  std::deque<double> a_p_memory = std::deque<double>(10, 0.0);

  // Initialize signals
  // set the driving cycle
  param->set_driving_cycle(data_v_h_106_p1, 0.1);

  size_t vector_size = size(param->get_driving_cycle());
  double loop_size = getDoubleValue(config["mpc"], "loop_size");
  if (loop_size < vector_size) {
    vector_size = loop_size;
  }

  std::vector<double> v_h_signal(vector_size, 0.0);
  std::vector<double> d_h_signal(vector_size, 0.0);
  std::vector<double> v_p_signal(vector_size, 0.0);
  std::vector<double> a_p_signal(vector_size, 0.0);
  std::vector<double> u_input_signal(vector_size, 0.0);
  std::vector<double> cost_level_signal(vector_size, 0.0);
  std::vector<double> num_iterations_signal(vector_size, 0.0);
  std::vector<double> cal_time_signal(vector_size, 0.0);
  std::vector<double> f_t_signal(vector_size, 0.0);
  std::vector<double> f_b_signal(vector_size, 0.0);
  std::vector<double> delta_s_signal(vector_size, 0.0);
  std::vector<double> a_p_prediction_vector(vector_size, 0.0);
  std::vector<double> run_times_signal(run_times, 0.0);

  int simulation_mode = getIntValue(
      config["mpc"],
      "simulation_mode"); // 1 one-time simulation, 2 closed loop simulation
  switch (simulation_mode) {
  case 1: // one-time simulation
  {
    oneTimeSimulationG2O(param, mpc, optimizer, numberOfIterations);
  } break;
  case 2: // closed loop simulation
  {
    int solvers_to_loop = 0;
    if (ampl_solver == 99) {
      solvers_to_loop = 8; // test all solvers
    }

    for (int counter = 0; counter <= solvers_to_loop; counter++) {

      double compute_time_ms = 0.01;
      double num_iterations = 0.0;
      double cost_level = 0.0;
      double u_input = 0.0;

      for (size_t k = 0; k < vector_size; k += 1) {

        update_controller_parameters(k, mpc, param, a_p_memory, u_input);
        if (true) {
          std::cout << "v_h: " << param->get_v_h_0()
                    << ", d_h: " << param->get_d_h_0()
                    << " , f_t_prev: " << param->get_f_t_prev()
                    << " , f_b_prev: " << param->get_f_b_prev()
                    << ", delta_s: " << param->get_delta_s() << " , k: " << k
                    << std::endl;

          // print the predicted v_p
          std::cout << "Predicted_v_p: ";
          for (size_t i = 0; i < param->get_vp_prediction().size(); i++) {
            std::cout << param->get_vp_prediction().at(i) << " ";
          }
          std::cout << std::endl;
        }

        // Extract the input after the optimization is done
        double num_iterations_memory = 0.0;

        for (int i = 0; i < run_times; i++) {
          u_input =
              g2o_compute_input(mpc, optimizer, numberOfIterations,
                                compute_time_ms, num_iterations, cost_level);
          run_times_signal[i] = compute_time_ms;
          if (i == 0) {
            num_iterations_memory = num_iterations;
          } else {
            /*       if (num_iterations != num_iterations_memory) {
                    // Throw an exception with a descriptive error message
                    throw std::runtime_error(
                        "Error: Number of iterations changed across runs. "
                        "Previous: " +
                        std::to_string(num_iterations_memory) +
                        ", Current: " + std::to_string(num_iterations));
                  } */
          }
        }

        v_h_signal[k] = param->get_v_h_0();
        d_h_signal[k] = param->get_d_h_0();
        //      std::cout << "d_h 0 " << param->get_d_h_0() << std::endl;
        v_p_signal[k] = param->get_vp_prediction()[0];
        a_p_signal[k] = a_p_memory.back();
        u_input_signal[k] = u_input;
        num_iterations_signal[k] = num_iterations;
        cal_time_signal[k] = meanWithoutOutliers(run_times_signal);
        cost_level_signal[k] = cost_level;
        f_t_signal[k] = param->get_f_t_prev();
        f_b_signal[k] = param->get_f_b_prev();
        delta_s_signal[k] = param->get_delta_s();
        a_p_prediction_vector[k] = compute_a_p_predicted(a_p_memory);

        // print_resutls(sol,var);

        std::cout << k << ": u_input = " << u_input << ", ";
        std::cout << "Number of iterations: " << num_iterations
                  << ", Computation time (ms): " << cal_time_signal[k] << " ms"
                  << std::endl;

        // Store results

        std::cout << "End of MPC "
                     "iteration, loop size = "
                  << k << " >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;

        if (num_iterations == 100) {
          //    break;
        }
      }
      // ################## save results and plot them ##################
      namespace fs = std::filesystem;
      // Save results to a file or plot them
      // create the file name:
      // linearSolverType_AlgEqType_mpcHorizon_"PIPM".bin
      std::time_t t = std::time(nullptr);
      std::tm tm = *std::localtime(&t);
      std::ostringstream time_stream;
      time_stream << std::put_time(&tm, "%y%m%d_%H%M");
      std::string time_string = time_stream.str();

      time_string += std::string("_G2O") + "_" + AlIneqType + "_" + AlgEqType +
                     "_" + solverType + "_" + std::to_string(linearSolverType);

      std::string file_name =
          time_string + "_" + std::to_string(mpcHorizon) + ".bin";
      std::cout << "File name: " << file_name << std::endl;

      // 1. Construct the target directory path using the CMake macro
      fs::path root = PROJECT_SOURCE_DIR;
      fs::path lib_dir =
          root / ".." / "lib"; // Assuming the library is in the "lib"
                               // subdirectory of the project root
      fs::path plotter_dir = root / "results";
      fs::path target_dir = plotter_dir / "data";

      // 2. Safety Check
      if (!fs::exists(target_dir)) {
        fs::create_directories(target_dir);
      }

      fs::path full_path = target_dir / file_name;

      DataSaver data_saver(full_path.string());

      saveSignals(data_saver, f_t_signal, f_b_signal, delta_s_signal,
                  a_p_prediction_vector);

      saveSignals(data_saver, v_h_signal, d_h_signal, v_p_signal, a_p_signal,
                  u_input_signal, num_iterations_signal, cal_time_signal,
                  cost_level_signal);

      std::string terminationCriterionStr;
      if (terminationCriterion == 0) {
        terminationCriterionStr = "UpdateNorm";
      } else if (terminationCriterion == 1) {
        terminationCriterionStr = "NewtonDecrement";
      } else if (terminationCriterion == 2) {
        terminationCriterionStr = "GradientNorm";
      }
      saveConfigSettings(data_saver, file_name, numberOfIterations, solverType,
                         AlgEqType, AlIneqType, terminationCriterionStr,
                         linearSolverType, mpcHorizon, ampl_solver);
      saveConfigSettings(data_saver, param, optimizer);

      data_saver.close();
      // run the python script to plot the results
      //  the python modular path

      fs::path script_path = plotter_dir / "results_harvesting.py";

      std::string command =
          "export PYTHONPATH=\"" + lib_dir.string() + ":$PYTHONPATH\"";
          
      command = command + " && python3 " + script_path.string() + " --dir " +
                target_dir.string() + " --data " + file_name;

      std::cout << command << std::endl;
      system(command.c_str());
      // plt::show();
    }
  }
  }

  return 0;
}
