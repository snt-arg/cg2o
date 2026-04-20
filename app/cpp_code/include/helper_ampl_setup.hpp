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

#ifndef HELPER_AMPL_SETUP_H_
#define HELPER_AMPL_SETUP_H_

#include "ampl/ampl.h"
#include "ampl/dataframe.h"
#include "ampl/entity.h"
#include <Eigen/Core>
#include <chrono>
#include <iostream>
#include <ranges>
#include <regex>
#include <string>

;
Eigen::VectorXd columnToEigenVector(const ampl::DataFrame &df,
                                    const char *colName) {

  /**
  Convert an AMPL DataFrame column to an Eigen vector.
  \param df The AMPL DataFrame containing the column.
  \param colName The name of the column to convert.
  \return An Eigen::VectorXd containing the values of the specified column.
  example usage:
  ampl::Parameter param = op_ampl.getParameter("v_p");
  auto out = columnToEigenVector(param.getValues(), param.name());

  */
  ampl::DataFrame::Column col = df.getColumn(colName);

  Eigen::VectorXd out(static_cast<Eigen::Index>(col.size()));
  for (Eigen::Index i = 0; i < out.size(); ++i) {
    out[i] = col[static_cast<std::size_t>(i)].dbl();
  }
  return out;
}

inline void printAmplVectorFlat(const ampl::Parameter &input) {
  std::cout << input.name() << ": "; // e.g. "v_p"

  const std::size_t n = static_cast<std::size_t>(input.numInstances());
  for (std::size_t i = 0; i < n - 1; ++i) {
    std::cout << input.get(i).dbl() << ", ";
  }
  std::cout << input.get(n - 1).dbl();
  std::cout << std::endl;
}

static double extract_iterations_from_solve_message(std::string &msg,
                                                    int ampl_solver = -1) {
  switch (ampl_solver) {
  case 0: { // ipopt
    std::regex pattern(R"(Number of Iterations\.+:\s*(\d+))");
    std::smatch matches;

    if (std::regex_search(msg, matches, pattern)) {
      if (matches.size() > 1) {
        return std::stod(matches[1].str());
      }
    }

    return -1.0; // not found
  }

  case 1:  // knitro with alg 1
  case 2:  // knitro   with alg 2
  case 3:  // knitro with alg 3
  case 4:  // knitro with alg 4
  case 23: // knitro with alg 5
  case 24: // knitro with alg 6
  {
    static const std::regex pat(R"(#\s*of\s*iterations\s*=\s*([+-]?\d+))",
                                std::regex::icase);

    std::smatch m;
    if (std::regex_search(msg, m, pat) && m.size() > 1) {
      return std::stod(m[1].str());
    }
    return -1.0;
    break;
  }
  case 5: // conopt
  {
    static const std::regex pat(R"(\b(\d+)\s+iterations?\b)",
                                std::regex::icase);

    std::smatch m;
    if (std::regex_search(msg, m, pat) && m.size() > 1) {
      return std::stod(m[1].str());
    }
    return -1.0;
    break;
  }
  case 6: // loqo
  {
    static const std::regex pat(R"(\(\s*(\d+)\s+iterations\b)",
                                std::regex::icase);

    std::smatch m;
    if (std::regex_search(msg, m, pat) && m.size() > 1) {
      return std::stod(m[1].str());
    }
    return -1.0;
    break;
  }
  case 7: // minos
  case 8: // snopt
  {
    static const std::regex pat(R"((\d+)\s+iterations\b)", std::regex::icase);

    std::smatch m;
    if (std::regex_search(msg, m, pat) && m.size() > 1) {
      return std::stod(m[1].str());
    }
    return -1.0;
    break;
  }
  case 9: { // bonmin
    // We capture the "It" field (25) which appears before the final time
    // field.
    std::regex r(R"(NLP0014I\s+\d+\s+\w+\s+[-+0-9.eE]+\s+(\d+)\s+[-+0-9.eE]+)");
    double iters = -1.0;

    for (std::sregex_iterator it(msg.begin(), msg.end(), r), end; it != end;
         ++it) {
      iters = std::stod((*it)[1].str()); // keep overwriting -> last match wins
    }
    return iters;
  } break;

  default:
    return -2; // not implemented
    break;
  }
}

void setSolverOptions(ampl::AMPL &op_ampl, int ampl_solver, int max_iter) {
  switch (ampl_solver) {
  case 0: // ipopt
    op_ampl.eval("option solver ipopt;");
    op_ampl.eval("option ipopt linear_solver HSL_MA77;");
    op_ampl.eval("option ipopt_options "
                 "'max_iter=" +
                 std::to_string(max_iter) + // maximum number of iterations
                 "';");                     // 1 ==> write .sol file

    // 2 ==> print primal variable values
    // 4 ==> print dual variable values
    // 8 ==> do not print solution message

    break;
  case 1:
    op_ampl.eval("option solver knitro;"); // KNITRO Optimization Solver
    op_ampl.eval("option knitro_options "
                 "'algorithm=1 " // 0 automatic, 1 Interior/Direct , 2
                                 // Interior/Conjugate-gradient 3 Active-set
                                 // algorithm, 4:   (SQP)
                 "maxit=" +
                 std::to_string(max_iter) + // maximum number of iterations
                 " numthreads=0"
                 " outlev=4';"); // 0 no output, 1 some output, 2 more output
    // simplex, 3 barrier
    break;
  case 2:
    op_ampl.eval("option solver knitro;"); // KNITRO Optimization Solver
    op_ampl.eval("option knitro_options "
                 "'algorithm=2 " // 0 automatic, 1 Interior/Direct , 2
                                 // Interior/Conjugate-gradient 3 Active-set
                                 // algorithm, 4:   (SQP)
                 "maxit=" +
                 std::to_string(max_iter) + // maximum number of iterations
                 " numthreads=0"
                 " outlev=4';"); // 0 no output, 1 some output, 2 more output
    // simplex, 3 barrier
    break;
  case 3:
    op_ampl.eval("option solver knitro;"); // KNITRO Optimization Solver
    op_ampl.eval("option knitro_options "
                 "'algorithm=3 " // 0 automatic, 1 Interior/Direct , 2
                                 // Interior/Conjugate-gradient 3 Active-set
                                 // algorithm, 4:   (SQP)
                 "maxit=" +
                 std::to_string(max_iter) + // maximum number of iterations
                 " numthreads=0"
                 " outlev=4';"); // 0 no output, 1 some output, 2 more output
    // simplex, 3 barrier
    break;
  case 4:
    op_ampl.eval("option solver knitro;"); // KNITRO Optimization Solver
    op_ampl.eval("option knitro_options "
                 "'algorithm=4 " // 0 automatic, 1 Interior/Direct , 2
                                 // Interior/Conjugate-gradient 3 Active-set
                                 // algorithm, 4:   (SQP)
                 "maxit=" +
                 std::to_string(max_iter) + // maximum number of iterations
                 " numthreads=0"
                 " outlev=4';"); // 0 no output, 1 some output, 2 more output
    // simplex, 3 barrier
    break;
  case 5:
    op_ampl.eval("option solver conopt;"); // CONOPT Nonlinear Optimization
    op_ampl.eval("option conopt_options 'threads =8"
                 " IterLim = 2500"
                 " outlev = 2';"); // 0 no output, 1 some output, 2 more output
    break;
  case 6:
    op_ampl.eval("option solver loqo;"); // LOQO Optimization Solver
    op_ampl.eval("option loqo_options 'maxit=" +
                 std::to_string(2 * max_iter) + // maximum number of iterations
                 "';"); // 0 no output, 1 some output, 2 more output
    break;
  case 7:
    op_ampl.eval("option solver minos;"); // MINOS Nonlinear Optimization
    op_ampl.eval("option minos_options 'iterations= 2500"
                 " feas_tol=1e-1"
                 "=';"); // 0 no output, 1 some output, 2 more output
    break;
  case 8:
    op_ampl.eval("option solver snopt;"); // SNOPT Optimization Solver
    op_ampl.eval("option snopt_options "
                 "' iterations=2500 feas_tol=1e-1"
                 "';");

    break;
  case 9:
    op_ampl.eval("option solver bonmin;");
    break;
  case 10:
    op_ampl.eval("option solver couenne;");
    break;
  case 11:
    op_ampl.eval("option solver highs;");
    break;
  case 12:
    op_ampl.eval("option solver scip;");
    break;
  case 13:
    op_ampl.eval("option solver gcg;");
    break;
  case 14:
    op_ampl.eval("option solver  cbc ;");
    break;
  case 15:
    op_ampl.eval("option solver gurobi;");
    // op_ampl.setOption("gurobi_options",
    //                 " Method = 2"); // Barrier method
    // Method = 0 → Primal Simplex  Method = 1
    // → Dual Simplex    Method = 2 →
    // Barrier             Method = -1 → Automatic(default)
    break;
  case 16:
    op_ampl.eval("option solver xpress;"); // Xpress Optimization Solver
    break;
  case 17:
    op_ampl.eval("option solver cplex;"); // CPLEX Optimization Solver
    break;
  case 18:
    op_ampl.eval("option solver copt;"); // COPT Optimization Solver
    break;
  case 19:
    op_ampl.eval("option solver mosek;"); // MOSEK Optimization Solver
    break;
  case 20:
    op_ampl.eval("option solver baron;"); // BARON Optimization Suite
    break;
  case 21:
    op_ampl.eval("option solver lgo;"); // LGO Optimization Solver
    break;
  case 22:
    op_ampl.eval(
        "option solver lindoglobal;"); // LINDO Global Optimization Solver
    break;
  default:
    std::runtime_error("Invalid solver selection.");
  }
}

void amplProblemSetup(ampl::AMPL &op_ampl, std::string ampl_files_path,
                      auto param, int ampl_solver, int max_iter) {

  ampl_files_path += "/mpc_ampl_files/";

  op_ampl.read(ampl_files_path + "/main.mod");
  op_ampl.read(ampl_files_path + "/cost.mod");
  op_ampl.read(ampl_files_path + "/eq_con.mod");
  op_ampl.read(ampl_files_path + "/ineq_con.mod");

  // setup the horizon
  ampl::Parameter N = op_ampl.getParameter("N");
  N.set(param->get_N());
  setSolverOptions(op_ampl, ampl_solver, max_iter);
}

double ampl_compute_input(auto &op_ampl, int ampl_solver, auto &amplP,
                          auto &param, double &compute_time_ms,
                          double &num_iterations, double &cost_level,
                          std::string &ampl_msg, int ampl_sovler,
                          int max_iter) {

  double v_h = param->get_v_h_0();
  double d_h = param->get_d_h_0();
  double f_t_prev = param->get_f_t_prev();
  double f_b_prev = param->get_f_b_prev();

  amplP.at("v_h_0").set(v_h);
  amplP.at("d_h_0").set(d_h);
  amplP.at("v_p").setValues(param->get_vp_prediction().data(),
                            param->get_vp_prediction().size());

  amplP.at("f_t_prev").set(f_t_prev);
  amplP.at("f_b_prev").set(f_b_prev);
  amplP.at("delta_s").set(param->get_delta_s());

  if (param->useSpaceDomain()) {
    amplP.at("space_mode").set(1.0); // space domain
  } else {
    amplP.at("space_mode").set(0.0); // time domain
  }

  op_ampl.setDblOption("presolve", 0);
  // op_ampl.setDblOption("max_time", 20.0); // set a time limit of 20
  // seconds

  setSolverOptions(op_ampl, ampl_solver, max_iter);
  op_ampl.eval("redeclare var v_h {K} := init_guess;");
  op_ampl.eval("redeclare var f_t {K_u} := init_guess;");
  op_ampl.eval("redeclare var f_b {K_u} := init_guess;");
  op_ampl.eval("redeclare var slack_1 {K_u} := init_guess;");
  op_ampl.eval("redeclare var slack_2 {K_u} := init_guess;");
  op_ampl.eval("redeclare var slack_3 {K_u} := init_guess;");

  auto start_time = std::chrono::high_resolution_clock::now();
  ampl_msg = op_ampl.getOutput("solve;");
  auto end_time = std::chrono::high_resolution_clock::now();
  compute_time_ms =
      std::chrono::duration<double, std::milli>(end_time - start_time).count();

  num_iterations = static_cast<int>(
      extract_iterations_from_solve_message(ampl_msg, ampl_solver));

  double u_input =
      op_ampl.getValue("f_t[0]").dbl() - op_ampl.getValue("f_b[0]").dbl();

  cost_level = op_ampl.getValue("cost").dbl();

  return u_input;
}

#endif // HELPER_AMPL_SETUP_H_