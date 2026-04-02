#ifndef MPC_FORMULATION_HPP
#define MPC_FORMULATION_HPP

#include "mpc_formulation.h"
#include "mpc_vertices.h"
#include <iostream>
#include <numeric>
#include <optimization_algorithm_property.h>
#include <ostream>
#include <stdexcept>

namespace cg2o::mpc {

template <typename Optimizer>
MPCFormulation<Optimizer>::MPCFormulation(int N,
                                          std::shared_ptr<Optimizer> optimizer,
                                          std::shared_ptr<MPCParameters> param)
    : _N(N), _optimizer(optimizer), _param(param) {
  computeOffsets();
}

template <typename Optimizer> void MPCFormulation<Optimizer>::computeOffsets() {
  _offset_d_h = 0 * (_N + 1);
  _offset_slack_1 = 1 * (_N + 1);
  _offset_f_t = 2 * (_N + 1);
  _offset_slack_2 = 3 * (_N + 1);
  _offset_v_h = 4 * (_N + 1);
  _offset_f_b = 5 * (_N + 1);
  _offset_slack_3 = 6 * (_N + 1);
  _offset_eqDynamics_v_h = 7 * (_N + 1);
  _offset_eqDynamics_d_h = 8 * (_N + 1);
}

template <typename Optimizer> void MPCFormulation<Optimizer>::setN(int N) {
  _N = N;
};

template <typename Optimizer> int MPCFormulation<Optimizer>::getN() const {
  return _N;
};

template <typename Optimizer>
void MPCFormulation<Optimizer>::setVLagrangianInitial(double vInitial) {
  _init_lagrange_multiplier = vInitial;
}

template <typename Optimizer>
double MPCFormulation<Optimizer>::getVLagrangianInitial() {
  return _init_lagrange_multiplier;
}

template <typename Optimizer>
void MPCFormulation<Optimizer>::addVertex(int offset, int k) {
  auto vertex = std::make_shared<Vertex_scalar>();
  vertex->setId(offset + k);
  // vertex->setId(offset + num_var_types * k); another way to set the id by
  // grouping the varaible in based on k
  //  if we have k continously increasingly num_var_types =7 and offset = 0,1,
  //  ... not 0*N, 1*N, ...
  _optimizer->addVertex(vertex.get());
  _vertices.push_back(vertex);
}
// add vertices

template <typename Optimizer>
void MPCFormulation<Optimizer>::addVertices(int k) {
  addVertex(_offset_v_h, k);
  addVertex(_offset_d_h, k);
  if (k == 0) {
    _optimizer->vertex(_offset_v_h + k)->setFixed(true);
    _optimizer->vertex(_offset_d_h + k)->setFixed(true);
  }

  if (k > 0) {
    addVertex(_offset_slack_1,
              k); // slack_1_0 does not exist because d_h_0 is fixed
  }
  if (k <= _N - 1) { // x_{N-1} is the last state
    addVertex(_offset_f_t, k);
    addVertex(_offset_f_b, k);
    addVertex(_offset_slack_2, k);
    addVertex(_offset_slack_3, k);
  }
}

template <typename Optimizer>
void MPCFormulation<Optimizer>::setVertexScalarEstimate(int id, double value) {
  // check if it is nan change it to safe value
  if (std::isnan(value)) {
    value = _replace_nan; // or any other safe default value
  }
  static_cast<Vertex_scalar *>(_optimizer->vertex(id))->setEstimate(value);
}

template <typename Optimizer>
double MPCFormulation<Optimizer>::getVertexScalarEstimate(int id) {
  try {
    auto vertex = _optimizer->vertex(id); // Direct call to optimizer
    if (!vertex) {
      throw std::runtime_error("Vertex with ID " + std::to_string(id) +
                               " not found.");
    }
    return static_cast<Vertex_scalar *>(vertex)->estimate();
  } catch (const std::runtime_error &e) {
    std::cerr << "Error in getVertexScalarEstimate: " << e.what() << std::endl;
    return std::numeric_limits<double>::quiet_NaN();
  }
}

// get the results
template <typename Optimizer>
std::vector<double> MPCFormulation<Optimizer>::getResults() {
  _results.clear();

  for (int k = 0; k <= _N; k++) {
    _results.push_back(getVertexScalarEstimate(_offset_v_h + k));
  }
  for (int k = 0; k <= _N; k++) {
    _results.push_back(getVertexScalarEstimate(_offset_d_h + k));
  }
  for (int k = 0; k <= _N - 1; k++) {
    _results.push_back(getVertexScalarEstimate(_offset_f_t + k));
  }
  for (int k = 0; k <= _N - 1; k++) {
    _results.push_back(getVertexScalarEstimate(_offset_f_b + k));
  }
  for (int k = 1; k <= _N; k++) {
    _results.push_back(getVertexScalarEstimate(_offset_slack_1 + k));
  }
  for (int k = 0; k <= _N - 1; k++) {
    _results.push_back(getVertexScalarEstimate(_offset_slack_2 + k));
  }
  for (int k = 0; k <= _N - 1; k++) {
    _results.push_back(getVertexScalarEstimate(_offset_slack_3 + k));
  }

  _results.push_back(0.0);
  _results.push_back(0.0);
  _results.push_back(0.0);

  for (int k = 0; k <= _N - 1; k++) {
    _results.push_back(getVertexScalarEstimate(_offset_eqDynamics_v_h + k));
  }

  _results.push_back(0.0);
  _results.push_back(0.0);
  _results.push_back(0.0);
  for (int k = 0; k <= _N - 1; k++) {
    _results.push_back(getVertexScalarEstimate(_offset_eqDynamics_d_h + k));
  }

  return _results;
}
// get the force input
template <typename Optimizer>
double MPCFormulation<Optimizer>::getForceInput(int k) {
  return getVertexScalarEstimate(_offset_f_t + k) -
         getVertexScalarEstimate(_offset_f_b + k);
}

template <typename Optimizer>
double MPCFormulation<Optimizer>::getVelocityDesired() {
  return getVertexScalarEstimate(_offset_v_h + 1);
}

template <typename Optimizer>
double MPCFormulation<Optimizer>::getAccelerationDesired(double v_h,
                                                         double road_grad) {
  // Evaluate the desired acceleration using the equation:
  // a = (u_car - m_total * g * sin(alpha) - f_roll - c_a * v_h^2) / m_eq
  auto u_car = this->getForceInput(0);
  // Get parameters
  auto m_total = _param->get_m_total();
  auto g = _param->get_g();       // Ensure g ≈ 9.81 m/s²
  auto c_a = _param->get_c_a();   // Aerodynamic drag coefficient
  auto m_eq = _param->get_m_eq(); // Equivalent mass (must be defined!)
  auto c_r = _param->get_c_r();   // Rolling resistance coefficient

  // Compute rolling resistance (always opposes motion)
  double f_roll;
  if (std::abs(v_h) > 0.05) {
    f_roll = c_r * m_total * g * cos(road_grad); // Note: Use alpha (radians)
  } else {
    f_roll = 0; // Or model static friction if needed
  }

  auto desired_acceleration =
      (u_car - m_total * g * sin(road_grad) - f_roll - c_a * v_h * v_h) / m_eq;

  return desired_acceleration;
}

// add edges
template <typename Optimizer>
void MPCFormulation<Optimizer>::addEdgeCost_brake_change(int k) {
  auto edge = std::make_shared<EdgeCost_brake_change>(k, _param);
  edge->setVertex(0, _optimizer->vertex(_offset_slack_3 + k));
  Eigen::Matrix<double, 1, 1> information;
  information << _param->get_BrakeChangePenalty();
  edge->setInformation(information);
  _optimizer->addEdge(edge.get());
  _edges.push_back(edge);
  _edges_cost.push_back(edge);
}

// add edges
template <typename Optimizer>
void MPCFormulation<Optimizer>::addEdgeCost_brake(int k) {
  auto edge = std::make_shared<EdgeCost_brake>(k, _param);
  edge->setVertex(0, _optimizer->vertex(_offset_f_b + k));
  Eigen::Matrix<double, 1, 1> information;
  information << _param->get_BrakePenalty();
  edge->setInformation(information);
  _optimizer->addEdge(edge.get());
  _edges.push_back(edge);
  _edges_cost.push_back(edge);
}

template <typename Optimizer>
void MPCFormulation<Optimizer>::addEdgeCost_energy(int k) {
  auto edge = std::make_shared<EdgeCost_energy>(k, _param);
  edge->setVertex(0, _optimizer->vertex(_offset_v_h + k));
  edge->setVertex(1, _optimizer->vertex(_offset_f_t + k));
  Eigen::Matrix<double, 2, 2> information;
  information << _param->get_Omega1() * edge->get_Q_c();
  edge->setInformation(information);
  _optimizer->addEdge(edge.get());
  _edges.push_back(edge);
  _edges_cost.push_back(edge);
}

template <typename Optimizer>
void MPCFormulation<Optimizer>::addEdgeCost_inner_distance(int k) {
  auto edge = std::make_shared<EdgeCost_inner_distance>(k, _param);
  edge->setVertex(0, _optimizer->vertex(_offset_slack_1 + k + 1));
  Eigen::Matrix<double, 1, 1> information;
  information << _param->get_Omega6();
  edge->setInformation(information);
  _optimizer->addEdge(edge.get());
  _edges.push_back(edge);
  _edges_cost.push_back(edge);
}

template <typename Optimizer>
void MPCFormulation<Optimizer>::addEdgeCost_traction_change(int k) {
  auto edge = std::make_shared<EdgeCost_traction_change>(k, _param);
  edge->setVertex(0, _optimizer->vertex(_offset_slack_2 + k));
  Eigen::Matrix<double, 1, 1> information;
  information << _param->get_TractionChangePenalty();
  edge->setInformation(information);
  _optimizer->addEdge(edge.get());
  _edges.push_back(edge);
  _edges_cost.push_back(edge);
}

template <typename Optimizer>
void MPCFormulation<Optimizer>::addEdgeDynamics_v_h(int k) {
  auto edge = std::make_shared<EdgeDynamics_v_h>(k, _param);
  edge->setVertex(0, _optimizer->vertex(_offset_v_h + k));
  edge->setVertex(1, _optimizer->vertex(_offset_v_h + k + 1));
  edge->setVertex(2, _optimizer->vertex(_offset_f_t + k));
  edge->setVertex(3, _optimizer->vertex(_offset_f_b + k));
  edge->setVertexLagrangeMultiplierId(_offset_eqDynamics_v_h + k);
  _optimizer->addEdgeEq(edge.get());
  _edges.push_back(edge);
}

template <typename Optimizer>
void MPCFormulation<Optimizer>::addEdgeDynamics_d_h(int k) {
  auto edge = std::make_shared<EdgeDynamics_d_h>(k, _param);
  edge->setVertex(0, _optimizer->vertex(_offset_d_h + k));
  edge->setVertex(1, _optimizer->vertex(_offset_d_h + k + 1));
  edge->setVertex(2, _optimizer->vertex(_offset_v_h + k));
  edge->setVertex(3, _optimizer->vertex(_offset_v_h + k + 1));
  edge->setVertexLagrangeMultiplierId(_offset_eqDynamics_d_h + k);
  _optimizer->addEdgeEq(edge.get());
  _edges.push_back(edge);
}

template <typename Optimizer>
void MPCFormulation<Optimizer>::addEdgeIneq_v_h_max(int k) {
  auto edge = std::make_shared<EdgeIneq_v_h_max>(k, _param);
  edge->setVertex(0, _optimizer->vertex(_offset_v_h + k + 1));
  _optimizer->addEdgeIneq(edge.get());
  _edges.push_back(edge);
}

template <typename Optimizer>
void MPCFormulation<Optimizer>::addEdgeIneq_v_h_min(int k) {
  auto edge = std::make_shared<EdgeIneq_v_h_min>(k, _param);
  edge->setVertex(0, _optimizer->vertex(_offset_v_h + k + 1));
  _optimizer->addEdgeIneq(edge.get());
  _edges.push_back(edge);
}

template <typename Optimizer>
void MPCFormulation<Optimizer>::addEdgeIneq_d_h_max(int k) {
  auto edge = std::make_shared<EdgeIneq_d_h_max>(k, _param);
  edge->setVertex(0, _optimizer->vertex(_offset_d_h + k + 1));
  edge->setVertex(1, _optimizer->vertex(_offset_v_h + k + 1));
  edge->setVertex(2, _optimizer->vertex(_offset_slack_1 + k + 1));
  _optimizer->addEdgeIneq(edge.get());
  _edges.push_back(edge);
}

template <typename Optimizer>
void MPCFormulation<Optimizer>::addEdgeIneq_d_h_min(int k) {
  auto edge = std::make_shared<EdgeIneq_d_h_min>(k, _param);
  edge->setVertex(0, _optimizer->vertex(_offset_d_h + k + 1));
  edge->setVertex(1, _optimizer->vertex(_offset_v_h + k + 1));
  _optimizer->addEdgeIneq(edge.get());
  _edges.push_back(edge);
}

template <typename Optimizer>
void MPCFormulation<Optimizer>::addEdgeIneq_f_b_change_limit(int k) {
  auto edge = std::make_shared<EdgeIneq_f_b_change_limit>(k, _param);
  if (k == 0) {
    edge->setVertex(0,
                    _optimizer->vertex(_offset_f_t +
                                       k)); // not needed in the first iteration
  } else {
    edge->setVertex(0, _optimizer->vertex(_offset_f_b + k - 1));
  }
  edge->setVertex(1, _optimizer->vertex(_offset_f_b + k));
  edge->setVertex(2, _optimizer->vertex(_offset_slack_3 + k));
  _optimizer->addEdgeIneq(edge.get());
  _edges.push_back(edge);
}

template <typename Optimizer>
void MPCFormulation<Optimizer>::addEdgeIneq_f_b_max(int k) {
  auto edge = std::make_shared<EdgeIneq_f_b_max>(k, _param);
  edge->setVertex(0, _optimizer->vertex(_offset_f_b + k));
  edge->setVertex(1, _optimizer->vertex(_offset_v_h + k));
  _optimizer->addEdgeIneq(edge.get());
  _edges.push_back(edge);
}

template <typename Optimizer>
void MPCFormulation<Optimizer>::addEdgeIneq_f_b_min(int k) {
  auto edge = std::make_shared<EdgeIneq_f_b_min>(k, _param);
  edge->setVertex(0, _optimizer->vertex(_offset_f_b + k));
  _optimizer->addEdgeIneq(edge.get());
  _edges.push_back(edge);
}

template <typename Optimizer>
void MPCFormulation<Optimizer>::addEdgeIneq_f_t_change_limit(int k) {
  auto edge = std::make_shared<EdgeIneq_f_t_change_limit>(k, _param);
  if (k == 0) {
    edge->setVertex(0,
                    _optimizer->vertex(_offset_f_b +
                                       k)); // not needed in the first iteration
                                            // so use any one of the vertices
  } else {
    edge->setVertex(0, _optimizer->vertex(_offset_f_t + k - 1));
  }
  edge->setVertex(1, _optimizer->vertex(_offset_f_t + k));
  edge->setVertex(2, _optimizer->vertex(_offset_slack_2 + k));
  _optimizer->addEdgeIneq(edge.get());
  _edges.push_back(edge);
}

template <typename Optimizer>
void MPCFormulation<Optimizer>::addEdgeIneq_f_t_max(int k) {
  auto edge = std::make_shared<EdgeIneq_f_t_max>(k, _param);
  edge->setVertex(0, _optimizer->vertex(_offset_f_t + k));
  edge->setVertex(1, _optimizer->vertex(_offset_v_h + k));
  _optimizer->addEdgeIneq(edge.get());
  _edges.push_back(edge);
}

template <typename Optimizer>
void MPCFormulation<Optimizer>::addEdgeIneq_f_t_min(int k) {
  auto edge = std::make_shared<EdgeIneq_f_t_min>(k, _param);
  edge->setVertex(0, _optimizer->vertex(_offset_f_t + k));
  _optimizer->addEdgeIneq(edge.get());
  _edges.push_back(edge);
}

template <typename Optimizer> void MPCFormulation<Optimizer>::addEdges(int k) {
  if (_optimizer->verbose()) {
    std::cout << "Adding edges for k = " << k << std::endl;
  }

  addEdgeCost_brake(k);        // f_b
  addEdgeCost_brake_change(k); // slack_3
  addEdgeCost_energy(k);
  addEdgeCost_inner_distance(k);  // min slack_1^2
  addEdgeCost_traction_change(k); // min slack_2^2

  #ifndef TEST_LINER_SOLVERS
  addEdgeDynamics_v_h(k);
  addEdgeDynamics_d_h(k);
  #endif

  addEdgeIneq_v_h_min(k);
  addEdgeIneq_v_h_max(k);

  addEdgeIneq_d_h_min(k);
  addEdgeIneq_d_h_max(k);

  addEdgeIneq_f_t_max(k);
  addEdgeIneq_f_t_min(k);

  addEdgeIneq_f_b_max(k);
  addEdgeIneq_f_b_min(k);

  addEdgeIneq_f_t_change_limit(k);
  addEdgeIneq_f_b_change_limit(k);
}

template <typename Optimizer>
bool MPCFormulation<Optimizer>::setInitialGuess(
    bool reset_lagrange_multipliers) {
  // set the initial guess for the optimization

  // uddate v_h_0 and d_h_0
  setVertexScalarEstimate(_offset_v_h, _param->get_v_h_0());
  setVertexScalarEstimate(_offset_d_h, _param->get_d_h_0());

  double shift = _param->get_interior_shift();
  double h_min = _param->get_h_min();
  double d_min = _param->get_d_min();
  double h_soft = _param->get_h_soft();
  double a2 = _param->get_a2();
  double b2 = _param->get_b2();
  double Brems_a1 = _param->get_Brems_a1();
  double Brems_b1 = _param->get_Brems_b1();
  double b_max = _param->get_b_max();
  double v_h_temp_km1;
  double f_t_k_prev;
  double f_b_k_prev;

  for (int k = 0; k <= _N - 1; k++) {
    if (k == 0) {
      v_h_temp_km1 = _param->get_v_h_0();
      f_t_k_prev = _param->get_f_t_prev();
      f_b_k_prev = _param->get_f_b_prev();
    } else {
      v_h_temp_km1 = getVertexScalarEstimate(_offset_v_h + k);
      f_t_k_prev = getVertexScalarEstimate(_offset_f_t + k - 1);
      f_b_k_prev = getVertexScalarEstimate(_offset_f_b + k - 1);
    }

    double v_h_temp_k = getVertexScalarEstimate(_offset_v_h + k + 1);
    v_h_temp_k =
        std::max(v_h_temp_k, _param->get_v_min_prediction(k + 1) + shift);
    v_h_temp_k =
        std::min(v_h_temp_k, _param->get_v_max_prediction(k + 1) - shift);
    setVertexScalarEstimate(_offset_v_h + k + 1, v_h_temp_k);

    double d_h_temp_k = getVertexScalarEstimate(_offset_d_h + k + 1);
    d_h_temp_k = std::max(d_h_temp_k, h_min * v_h_temp_k + d_min + shift);
    setVertexScalarEstimate(_offset_d_h + k + 1, d_h_temp_k);

    double slack_1_temp_k = getVertexScalarEstimate(_offset_slack_1 + k + 1);
    slack_1_temp_k =
        std::max(slack_1_temp_k, d_h_temp_k - h_soft * v_h_temp_k + shift);
    setVertexScalarEstimate(_offset_slack_1 + k + 1, slack_1_temp_k);

    double f_t_temp_k = getVertexScalarEstimate(_offset_f_t + k);
    f_t_temp_k = std::max(f_t_temp_k, 0 + shift);
    f_t_temp_k = std::min(f_t_temp_k, a2 * v_h_temp_km1 + b2 - shift);
    setVertexScalarEstimate(_offset_f_t + k, f_t_temp_k);

    double slack_2_temp_k = getVertexScalarEstimate(_offset_slack_2 + k);
    slack_2_temp_k = std::max(slack_2_temp_k, f_t_temp_k - f_t_k_prev + shift);
    slack_2_temp_k =
        std::max(slack_2_temp_k, -(f_t_temp_k - f_t_k_prev) + shift);
    setVertexScalarEstimate(_offset_slack_2 + k, slack_2_temp_k);

    double f_b_temp_k = getVertexScalarEstimate(_offset_f_b + k);
    f_b_temp_k = std::max(f_b_temp_k, 0 + shift);
    f_b_temp_k = std::min(f_b_temp_k,
                          Brems_a1 * v_h_temp_km1 + (b_max + Brems_b1) - shift);
    setVertexScalarEstimate(_offset_f_b + k, f_b_temp_k);

    double slack_3_temp_k = getVertexScalarEstimate(_offset_slack_3 + k);
    slack_3_temp_k = std::max(slack_3_temp_k, f_b_temp_k - f_b_k_prev + shift);
    slack_3_temp_k =
        std::max(slack_3_temp_k, -(f_b_temp_k - f_b_k_prev) + shift);
    setVertexScalarEstimate(_offset_slack_3 + k, slack_3_temp_k);
  }

  if (reset_lagrange_multipliers) {
    _optimizer->resetLagrangeMultiplierEq(); // automatically done by the
                                             // optimizer but we do it to see
                                             // the initial values
  }
  return true;
}

template <typename Optimizer>
bool MPCFormulation<Optimizer>::setFixedInitialGuess(double value) {
  // set the initial guess for the optimization

  // uddate v_h_0 and d_h_0
  setVertexScalarEstimate(_offset_v_h, _param->get_v_h_0());
  setVertexScalarEstimate(_offset_d_h, _param->get_d_h_0());

  for (int k = 0; k <= _N - 1; k++) {
    setVertexScalarEstimate(_offset_v_h + k + 1, value);
    setVertexScalarEstimate(_offset_d_h + k + 1, value);
    setVertexScalarEstimate(_offset_f_t + k, value);
    setVertexScalarEstimate(_offset_f_b + k, value);
    setVertexScalarEstimate(_offset_slack_1 + k + 1, value);
    setVertexScalarEstimate(_offset_slack_2 + k, value);
    setVertexScalarEstimate(_offset_slack_3 + k, value);
  }

  return true;
}

template <typename Optimizer> void MPCFormulation<Optimizer>::setupMPC() {
  std::cout << "Setting up MPC formulation with Horizon " << _N << std::endl;

  // create the prediction vector for alpha, v_min, v_man (can be done later
  // in the predictor)
  _param->set_alpha_prediction(_N + 1, _param->get_alpha());
  _param->set_v_min_prediction(_N + 1, _param->get_v_min());
  _param->set_v_max_prediction(_N + 1, _param->get_v_max());

  // creeate verteices
  for (int k = 0; k <= _N; k++) {
    addVertices(k);
  }

  // create edges

  for (int k = 0; k <= _N - 1; k++) {
    addEdges(k);
  }
}

// For closed loop simulation
template <typename Optimizer>
double MPCFormulation<Optimizer>::simulatorVelocity(double u_car) {
  if (u_car < -7000) {
    u_car = -7000; // ensure the pyhsical limits of the car
  }
  double v_h_0 = _param->get_v_h_0();
  _param->set_v_h_prev(v_h_0); // update the previous velocity

  double alpha =
      _param->get_alpha_prediction(0); // Gradient in radians, assumed flat
  double Ts = _param->get_Ts();
  double g = _param->get_g();
  double m_total = _param->get_m_v() + _param->get_m_p();
  double m_eq = m_total + _param->get_e_i() * _param->get_m_v();
  double c_r = _param->get_x1() * _param->get_c_r();
  double c_a = _param->get_x2() * _param->get_c_a();

  double f_roll = (v_h_0 > 0.05) ? (c_r * m_total * g * cos(alpha)) : 0.0;

  // Calculate new velocity
  double v_h_1 = (Ts / m_eq) * (u_car - m_total * g * sin(alpha) - f_roll -
                                c_a * v_h_0 * v_h_0) +
                 v_h_0;
  v_h_1 = std::max(v_h_1, 0.0); // Ensure velocity is non-negative
  _param->set_v_h_0(v_h_1); // update the initial velocity

  return v_h_1;
}

// For closed loop simulation
template <typename Optimizer>
void MPCFormulation<Optimizer>::sensor(const double &v_p,
                                       std::deque<double> &a_p_memory) {
  // Measure the distance to preceding and its velocity
  double Ts = _param->get_Ts();

  double v_h = _param->get_v_h_0();
  double v_h_prev = _param->get_v_h_prev();

  double d_h = _param->get_d_h_0();
  double v_p_prev = _param->get_v_p_prev();

  d_h = d_h + Ts * (v_p_prev + v_p) / 2 - Ts * (v_h_prev + v_h) / 2;
  double a_p = (v_p - v_p_prev) / Ts;
  a_p_memory.push_back(a_p);
  a_p_memory.pop_front();

  // update the previous values
  _param->set_d_h_0(d_h);
  _param->set_v_p_prev(v_p);
}

template <typename Optimizer>
std::vector<double>
MPCFormulation<Optimizer>::predictor(double v_p,
                                     const std::deque<double> &a_p_memory) {
  std::vector<double> weight_accel_memory = {0, 0, 0, 0,    0,
                                             0, 0, 0, 0.25, 0.75};
  double predicted_accel =
      std::inner_product(weight_accel_memory.begin(), weight_accel_memory.end(),
                         a_p_memory.begin(), 0.0);

  double a_p = a_p_memory.back();
  std::vector<double> F(_N + 1);
  std::iota(F.begin(), F.end(), 0); // F = [0, 1, 2, ..., N]
  std::vector<double> v_p_prediction(F.size(), 0);
  if (predicted_accel >= 0) {
    predicted_accel = std::min(a_p, 0.5 * predicted_accel);
    for (size_t i = 0; i < F.size(); ++i) {
      v_p_prediction[i] = predicted_accel * F[i] + v_p;
      v_p_prediction[i] = predicted_accel * F[i] + v_p;
      if (v_p_prediction[i] < 0)
        v_p_prediction[i] = 0;

      if (i >= 3)
        v_p_prediction[i] = v_p_prediction[2]; // Fix values after the 3rd entry
    }
  } else {
    predicted_accel = std::max({a_p, predicted_accel, -0.2});
    for (size_t i = 0; i < F.size(); ++i) {
      v_p_prediction[i] = predicted_accel * F[i] + v_p;
      if (v_p_prediction[i] < 0)
        v_p_prediction[i] = 0;
    }
  }
  return v_p_prediction;
}

template <typename Optimizer> double MPCFormulation<Optimizer>::computeCost() {
  double cost = 0.0;

  for (const auto &edge : _edges_cost) {
    // Accumulate chi-squared error
    cost += edge->chi2(); // Assuming chi2() returns the squared error for the
                          // edge
  }

  return cost;
}

} // namespace cg2o::mpc
#endif // MPC_FORMULATION_HPP
