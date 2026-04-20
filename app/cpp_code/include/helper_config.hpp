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

#ifndef HELPER_CONFIG_H_
#define HELPER_CONFIG_H_

#include <filesystem>
#include <yaml-cpp/yaml.h>

static std::string getConfigPath(int argc, char **argv) {
  // 1) If user placed config.yaml in the run directory, use it
  if (std::filesystem::exists("config.yaml")) {
    return "config.yaml";
  }

  // 2) Allow optional override path
  for (int i = 1; i + 1 < argc; ++i) {
    if (std::string(argv[i]) == "--config") {
      return argv[i + 1];
    }
  }

  // 3) Fallback to CMake default
#ifdef DEFAULT_CONFIG_PATH
  return std::string(DEFAULT_CONFIG_PATH);
#else
  return "config/config.yaml";
#endif
}
inline int getIntValue(const YAML::Node &node, const std::string &key) {
  if (node[key]) {
    return node[key].as<int>();
  }
  throw std::runtime_error("Missing required YAML key: '" + key + "'");

  return 0;
}

inline double getDoubleValue(const YAML::Node &node, const std::string &key) {
  if (node[key]) {
    return node[key].as<double>();
  }
  throw std::runtime_error("Missing required YAML key: '" + key + "'");
  return 0;
}

inline std::string getStringValue(const YAML::Node &node,
                                  const std::string &key) {
  if (node[key]) {
    return node[key].as<std::string>();
  }
  throw std::runtime_error("Missing required YAML key: '" + key + "'");
  return "";
}

void configureMPCParameters(const YAML::Node &config, auto &param,
                            auto &optimizer) {

  int spaceDomain = getIntValue(config["acc_settings"], "spaceDomain");
  int linearDynamics = getIntValue(config["acc_settings"], "linearDynamics");
  int linearInequalities =
      getIntValue(config["acc_settings"], "linearInequalities");
  int scaleInequalities =
      getIntValue(config["acc_settings"], "scaleInequalities");
  param->set_space_domain(spaceDomain);
  param->set_linearDynamics(linearDynamics);
  param->set_linearInequalities(linearInequalities);
  param->set_scaleInequalities(scaleInequalities);
  // Initialize optimizer

  optimizer->setVerbose(getIntValue(config["general_settings"], "verbose"));
  optimizer->setNumInnerIterationsMax(
      getIntValue(config["general_settings"], "num_inner_iterations_max"));
  double alpha_backtracking =
      getDoubleValue(config["general_settings"], "alpha_backtracking");
  optimizer->setAlphaBacktracking(alpha_backtracking);
#ifdef USE_AL
  optimizer->setAlgorithmEq(AlgEqType);
  optimizer->_lagrange_multiplier_initial_ineq = getDoubleValue(
      config["general_settings"], "lagrange_multiplier_initial_eq");
  optimizer->_rho_initial_ineq =
      getDoubleValue(config["general_settings"], "rho_initial_eq");
  optimizer->_rho_min_ineq =
      getDoubleValue(config["general_settings"], "rho_min_eq");
  optimizer->_rho_max_ineq =
      getDoubleValue(config["general_settings"], "rho_max_eq");
  optimizer->_rho_update_factor_ineq =
      getDoubleValue(config["general_settings"], "rho_update_factor_eq");

#endif
#ifdef USE_BIPM
  optimizer->setKappaInitial(
      getDoubleValue(config["bipm_settings"], "kappa_initial"));
  optimizer->setKappaFinal(
      getDoubleValue(config["bipm_settings"], "kappa_final"));
  optimizer->setKappaUpdateFactor(
      getDoubleValue(config["bipm_settings"], "kappa_update_factor"));
#endif
#ifdef USE_ISPD
   optimizer->setStepSizeStrategy(
      getIntValue(config["ispd_settings"], "step_size_strategy"));
  optimizer->setIneqBacktrackingStepMin(
      getDoubleValue(config["ispd_settings"], "ineq_backtracking_step_min"));
  optimizer->setAuxBacktrackingStepMin(
      getDoubleValue(config["ispd_settings"], "aux_backtracking_step_min"));

  optimizer->setInitKappaStrategy(
      getIntValue(config["ispd_settings"], "_init_kappa_strategy"));
  optimizer->setKappaInitial(
      getDoubleValue(config["ispd_settings"], "kappa_initial"));

  optimizer->setUpdateKappaStrategy(
      getDoubleValue(config["ispd_settings"], "update_kappa_strategy"));
  optimizer->setKappaUpdateFactor(
      getDoubleValue(config["ispd_settings"], "kappa_update_factor"));
  optimizer->setKappaFinal(
      getDoubleValue(config["ispd_settings"], "kappa_final"));
  optimizer->setTau(
      getDoubleValue(config["ispd_settings"], "tau"));
  optimizer->setLimitKappaFinal(
      getIntValue(config["ispd_settings"], "limit_kappa_final"));

  optimizer->setInitSlackStrategy(
      getIntValue(config["ispd_settings"], "init_slack_strategy"));
  optimizer->setSlackVariableInitialIneq(
      getDoubleValue(config["ispd_settings"], "slack_variable_initial_ineq"));

  optimizer->setInitLagrangeStrategy(
      getIntValue(config["ispd_settings"], "init_lagrange_strategy"));
  optimizer->setLagrangeMultiplierInitialIneq(
      getDoubleValue(config["ispd_settings"], "lagrange_multiplier_initial_ineq"));

  optimizer->setIneqPredictionParam(
      getDoubleValue(config["ispd_settings"], "ineq_prediction_param"));
  optimizer->setIneqPredictionStepStrategy(
      getIntValue(config["ispd_settings"], "ineq_prediction_step_strategy"));
  optimizer->setIneqPredictionStrategy(
      getIntValue(config["ispd_settings"], "ineq_prediction_strategy"));

  optimizer->setKeepAuxPositiveStrategy(
      getIntValue(config["ispd_settings"], "keep_aux_positive_strategy"));
  optimizer->setAuxScalingFactor(
      getDoubleValue(config["ispd_settings"], "aux_scaling_factor"));
  optimizer->setAuxCorrectionValue(
      getDoubleValue(config["ispd_settings"], "aux_correction_value"));

#endif
}

void saveConfigSettings(auto &data_saver, auto &param, const auto &optimizer) {
  data_saver.write("alpha", optimizer->alpha());
  data_saver.write("IsSpaceDomain", param->isSpaceDomain());
  data_saver.write("useLinearDynamics", param->isLinearDynamics());
  data_saver.write("useLinearInequalities", param->isLinearInequalities());
  data_saver.write("scaleInequalities", param->isScaleInequalities());

#if defined(USE_AL)
  data_saver.write("init_lagrange_multiplier_ineq",
                   optimizer->_lagrange_multiplier_initial_ineq);
  data_saver.write("init_rho_ineq", optimizer->_rho_initial_ineq);
  data_saver.write("rho_min_ineq", optimizer->_rho_min_ineq);
  data_saver.write("rho_max_ineq", optimizer->_rho_max_ineq);
  data_saver.write("rho_upate_factor_ineq", optimizer->_rho_update_factor_ineq);
#endif
#ifdef USE_BIPM
  data_saver.write("kappa_final", optimizer->kappaFinal());
  data_saver.write("t_init", optimizer->kappaInitial());
  data_saver.write("nu_update", optimizer->kappaUpdateFactor());
#endif
#ifdef USE_ISPD
  data_saver.write("t_init", optimizer->kappaInitial());
  data_saver.write("nu_update", optimizer->kappaUpdateFactor());
  data_saver.write("step_size_strategy", optimizer->stepSizeStrategy());
  data_saver.write("ineq_backtracking_step_min",
                   optimizer->ineqBacktrackingStepMin());
  data_saver.write("aux_backtracking_step_min",
                   optimizer->auxBacktrackingStepMin());
  data_saver.write("_init_kappa_strategy", optimizer->initKappaStrategy());
  data_saver.write("kappa_initial", optimizer->kappaInitial());
  data_saver.write("update_kappa_strategy", optimizer->updateKappaStrategy());
  data_saver.write("nu", optimizer->kappaUpdateFactor());
  data_saver.write("kappa_final", optimizer->kappaFinal());
  data_saver.write("tau", optimizer->tau());
  data_saver.write("limit_kappa_final", optimizer->limitKappaFinal());
  data_saver.write("init_slack_strategy", optimizer->initSlackStrategy());
  data_saver.write("slack_variable_initial_ineq",
                   optimizer->slackVariableInitialIneq());
  data_saver.write("lagrange_multiplier_initial_ineq",
                   optimizer->lagrangeMultiplierInitialIneq());
  data_saver.write("ineq_prediction_strategy",
                   optimizer->ineqPredictionStrategy());
  data_saver.write("ineq_prediction_step_strategy",
                   optimizer->ineqPredictionStepStrategy());
  data_saver.write("ineq_prediction_param", optimizer->ineqPredictionParam());
  data_saver.write("keep_aux_positive_strategy",
                   optimizer->keepAuxPositiveStrategy());
  data_saver.write("aux_scaling_factor", optimizer->auxScalingFactor());
  data_saver.write("aux_correction_value", optimizer->auxCorrectionValue());
#endif
}

void saveSignals(auto &data_saver, std::vector<double> &v_h_signal,
                 std::vector<double> &d_h_signal,
                 std::vector<double> &v_p_signal,
                 std::vector<double> &a_p_signal,
                 std::vector<double> &u_input_signal,
                 std::vector<double> &num_iterations_signal,
                 std::vector<double> &cal_time_signal,
                 std::vector<double> &cost_level_signal) {

  data_saver.write("v_h_signal", v_h_signal);
  data_saver.write("d_h_signal", d_h_signal);
  data_saver.write("v_p_signal", v_p_signal);
  data_saver.write("a_p_signal", a_p_signal);
  data_saver.write("u_input_signal", u_input_signal);
  data_saver.write("num_iterations_signal", num_iterations_signal);
  data_saver.write("cal_time_signal", cal_time_signal);
  data_saver.write("cost_level_signal", cost_level_signal);
};

void saveSignals(auto &data_saver, std::vector<double> &f_t_signal,
                 std::vector<double> &f_b_signal,
                 std::vector<double> &delta_s_signal,
                 std::vector<double> &a_p_prediction_vector) {
  data_saver.write("f_t_signal", f_t_signal);
  data_saver.write("f_b_signal", f_b_signal);
  data_saver.write("delta_s_signal", delta_s_signal);
  data_saver.write("a_p_prediction_vector", a_p_prediction_vector);
};

void saveSignals(auto &data_saver, std::vector<std::string> &ampl_msg_signal) {
  data_saver.write("ampl_msg_signal", ampl_msg_signal);
};

void saveSignalsAmpl(auto &data_saver, std::vector<double> &u_input_signal_lang,
                     std::vector<double> &num_iterations_signal_lang,
                     std::vector<double> &cal_time_signal_lang,
                     std::vector<double> &cost_level_signal_lang) {
  data_saver.write("u_input_signal_lang", u_input_signal_lang);
  data_saver.write("num_iterations_signal_lang", num_iterations_signal_lang);
  data_saver.write("cal_time_signal_lang", cal_time_signal_lang);
  data_saver.write("cost_level_signal_lang", cost_level_signal_lang);
};

void saveConfigSettings(auto &data_saver, std::string &file_name,
                        int numberOfIterations, std::string &solverType,
                        std::string &AlgEqType, std::string &AlIneqType,
                        std::string &terminationCriterionStr,
                        int linearSolverType, int mpcHorizon, int ampl_solver) {
  data_saver.write("File_name", file_name);
  data_saver.write("iteration_number_max", numberOfIterations);
  data_saver.write("Alg_unconstrainted", solverType);
  data_saver.write("Alg_equality", AlgEqType);
  data_saver.write("Alg_inequality", AlIneqType);
  data_saver.write("Termination_criteria:", terminationCriterionStr);
  data_saver.write("linearSolverType", linearSolverType);
  data_saver.write("mpcHorizon", mpcHorizon);
  data_saver.write("ampl_solver", ampl_solver);
}

#endif // HELPER_CONFIG_H_