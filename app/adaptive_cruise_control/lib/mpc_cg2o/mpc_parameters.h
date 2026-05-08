#ifndef MPC_PARAMETERS_H
#define MPC_PARAMETERS_H

#include "driving_cycles.h"
#include <vector>

namespace cg2o::mpc {

class MPCParameters {
private:
  // Parameters as const private member
  const bool space_domain = true;
  const bool linearDynamics = false;
  const bool linearInequalities = false;
  const bool scaleInequalities = false;
  const bool use_space_domain = true;
  const double domain_switch_threshold =
    3.0; // m/s threshold between time and space domain
  const int N = 2;
  const int N_max = 10;
  const double interior_shift = 0.01;
  const double c_a = 0.421022;
  const double c_r = 0.00708572690577572;
  const double e_i = 0.07632073023674692;
  const double g = 9.81;
  const double m_eq = 1639.0408163265306;
  const double m_p = 200.0;
  const double m_total = 1537.0;
  const double m_v = 1337.0;
  const double p1 = 0.0;
  const double p2 = 20.0;
  const double x1 = 1.1;
  const double x2 = 2.0;
  const double Omega1 = 0.01;            // energy consumption weight
  const double Omega6 = 1000.0;          // distance weight
  const double BrakeChangePenalty = 0.1; // penalty for brake change
  const double BrakePenalty = 0.08;
  const double TractionChangePenalty = 0.01;
  const double Brems_a1 = -47.01763813449595;
  const double Brems_b1 = 3068.0;
  const double Regen_max = 4000.0;
  const double a1 = 0.0;
  const double a2 = -206.46604377690898;
  const double b1 = 0.0;
  const double b2 = 9669.0;
  const double b3 = 0.0;
  const double b4 = 6903.0;
  const double b_max = 5000.0;
  const double d_min = 4.5;
  const double delta_t = 1.0;
  const double delta_s = 1.00;
  const double h_min = 0.75;
  const double h_soft = 1.25;
  const double delta_f_t_max = 0;
  const double delta_f_b_max = 0;
  const std::vector<double> eff_map = {
      67.828402111868030,
      0.768901253903272,
      2.453536693183027e-05,
      -0.767327621762776,
      3.831640131729236e-04,
      0.001505129345850}; // p00= eff_map[0]; p01 = eff_map[1]; p02 =
                          // eff_map[2]; double p10 = eff_map[3]; p11 =
                          // eff_map[4]; p20 = eff_map[5];

  const double Ts = 0.1;
  const std::vector<double> driving_cycle = data_v_h_106_p1;
  const std::vector<double> vp_prediction = {0, 0};
  const std::vector<double> alpha_prediction = {0, 0};
  const std::vector<double> v_min_prediction = {0, 0};
  const std::vector<double> v_max_prediction = {0, 0};

  const double alpha = 0.0;
  const double v_min = -1.0;
  const double v_max = 35.0;
  const double f_t_prev = 0;
  const double f_b_prev = 0;
  const double v_p_prev = 0;
  const double v_h_prev = 0;
  const double v_h_0 = 0.0;
  const double d_h_0 = 4.5;

public:
  // Getters for read-only access
  bool isLinearDynamics() const { return linearDynamics; }
  bool isSpaceDomain() const { return space_domain; }
  bool isLinearInequalities() const { return linearInequalities; }
  bool isScaleInequalities() const { return scaleInequalities; }
  bool useSpaceDomain() const { return use_space_domain; }
  double get_domain_switch_threshold() const { return domain_switch_threshold; }
  int get_N() const { return N; }
  int get_N_max() const { return N_max; }
  double get_interior_shift() const { return interior_shift; }
  double get_Ts() const { return Ts; }
  double get_c_a() const { return c_a; }
  double get_c_r() const { return c_r; }
  double get_e_i() const { return e_i; }
  double get_g() const { return g; }
  double get_m_eq() const { return m_eq; }
  double get_m_p() const { return m_p; }
  double get_m_total() const { return m_total; }
  double get_m_v() const { return m_v; }
  double get_p1() const { return p1; }
  double get_p2() const { return p2; }
  double get_x1() const { return x1; }
  double get_x2() const { return x2; }
  double get_BrakeChangePenalty() const { return BrakeChangePenalty; }
  double get_BrakePenalty() const { return BrakePenalty; }
  double get_Brems_a1() const { return Brems_a1; }
  double get_Brems_b1() const { return Brems_b1; }
  double get_Omega1() const { return Omega1; }
  double get_Omega6() const { return Omega6; }
  double get_Regen_max() const { return Regen_max; }
  double get_TractionChangePenalty() const { return TractionChangePenalty; }
  double get_a1() const { return a1; }
  double get_a2() const { return a2; }
  double get_b1() const { return b1; }
  double get_b2() const { return b2; }
  double get_b3() const { return b3; }
  double get_b4() const { return b4; }
  double get_b_max() const { return b_max; }
  double get_d_min() const { return d_min; }
  double get_delta_t() const { return delta_t; }
  double get_delta_s() const { return delta_s; }
  double get_h_min() const { return h_min; }
  double get_h_soft() const { return h_soft; }
  double get_delta_f_t_max() const { return delta_f_t_max; }
  double get_delta_f_b_max() const { return delta_f_b_max; }
  std::vector<double> get_eff_map() const { return eff_map; }
  double get_driving_cycle(int k) const { return driving_cycle[k]; }
  double get_vp_prediction(int k) const { return vp_prediction[k]; }
  double get_alpha_prediction(int k) const { return alpha_prediction[k]; }
  double get_v_min_prediction(int k) const { return v_min_prediction[k]; }
  double get_v_max_prediction(int k) const { return v_max_prediction[k]; }
  std::vector<double> get_driving_cycle() const { return driving_cycle; }
  std::vector<double> get_vp_prediction() const { return vp_prediction; }
  std::vector<double> get_alpha_prediction() const { return alpha_prediction; }
  std::vector<double> get_v_min_prediction() const { return v_min_prediction; }
  std::vector<double> get_v_max_prediction() const { return v_max_prediction; }
  double get_alpha() const { return alpha; }
  double get_v_min() const { return v_min; }
  double get_v_max() const { return v_max; }
  double get_f_t_prev() const { return f_t_prev; }
  double get_f_b_prev() const { return f_b_prev; }
  double get_v_p_prev() const { return v_p_prev; }
  double get_v_h_prev() const { return v_h_prev; }
  double get_v_h_0() const { return v_h_0; }
  double get_d_h_0() const { return d_h_0; }

  // Setters for controlled modification (using const_cast)
  void set_linearDynamics(bool value) {
    const_cast<bool &>(linearDynamics) = value;
  }
  void set_space_domain(bool value) {
    const_cast<bool &>(space_domain) = value;
  }
  void set_linearInequalities(bool value) {
    const_cast<bool &>(linearInequalities) = value;
  }
  void set_scaleInequalities(bool value) {
    const_cast<bool &>(scaleInequalities) = value;
  }
  void set_use_space_domain(bool value) {
    const_cast<bool &>(use_space_domain) = value;
  }
  void set_domain_switch_threshold(double value) {
    const_cast<double &>(domain_switch_threshold) = value;
  }
  void set_N(int value) { const_cast<int &>(N) = value; }
  void set_N_max(int value) { const_cast<int &>(N_max) = value; }
  void set_interior_shift(double value) {
    const_cast<double &>(interior_shift) = value;
  }
  void set_Ts(double value) { const_cast<double &>(Ts) = value; }
  void set_c_a(double value) { const_cast<double &>(c_a) = value; }
  void set_c_r(double value) { const_cast<double &>(c_r) = value; }
  void set_e_i(double value) { const_cast<double &>(e_i) = value; }
  void set_m_eq(double value) { const_cast<double &>(m_eq) = value; }
  void set_m_p(double value) { const_cast<double &>(m_p) = value; }
  void set_m_total(double value) { const_cast<double &>(m_total) = value; }
  void set_m_v(double value) { const_cast<double &>(m_v) = value; }
  void set_p1(double value) { const_cast<double &>(p1) = value; }
  void set_p2(double value) { const_cast<double &>(p2) = value; }
  void set_x1(double value) { const_cast<double &>(x1) = value; }
  void set_x2(double value) { const_cast<double &>(x2) = value; }
  void set_BrakeChangePenalty(double value) {
    const_cast<double &>(BrakeChangePenalty) = value;
  }
  void set_BrakePenalty(double value) {
    const_cast<double &>(BrakePenalty) = value;
  }
  void set_Brems_a1(double value) { const_cast<double &>(Brems_a1) = value; }
  void set_Brems_b1(double value) { const_cast<double &>(Brems_b1) = value; }
  void set_Omega1(double value) { const_cast<double &>(Omega1) = value; }
  void set_Omega6(double value) { const_cast<double &>(Omega6) = value; }
  void set_Regen_max(double value) { const_cast<double &>(Regen_max) = value; }
  void set_TractionChangePenalty(double value) {
    const_cast<double &>(TractionChangePenalty) = value;
  }
  void set_a1(double value) { const_cast<double &>(a1) = value; }
  void set_a2(double value) { const_cast<double &>(a2) = value; }
  void set_b1(double value) { const_cast<double &>(b1) = value; }
  void set_b2(double value) { const_cast<double &>(b2) = value; }
  void set_b3(double value) { const_cast<double &>(b3) = value; }
  void set_b4(double value) { const_cast<double &>(b4) = value; }
  void set_b_max(double value) { const_cast<double &>(b_max) = value; }
  void set_d_min(double value) { const_cast<double &>(d_min) = value; }
  void set_delta_t(double value) { const_cast<double &>(delta_t) = value; }
  void set_delta_s(double value) { const_cast<double &>(delta_s) = value; }
  void set_h_min(double value) { const_cast<double &>(h_min) = value; }
  void set_h_soft(double value) { const_cast<double &>(h_soft) = value; }
  void set_delta_f_t_max(double value) {
    const_cast<double &>(delta_f_t_max) = value;
  }
  void set_delta_f_b_max(double value) {
    const_cast<double &>(delta_f_b_max) = value;
  }
  void set_eff_map(std::vector<double> value) {
    const_cast<std::vector<double> &>(eff_map) = value;
  }

  void set_driving_cycle(std::size_t size, double value) {
    const_cast<std::vector<double> &>(driving_cycle).assign(size, value);
  }
  void set_alpha_prediction(std::size_t size, double value) {
    const_cast<std::vector<double> &>(alpha_prediction).assign(size, value);
  }
  void set_v_min_prediction(std::size_t size, double value) {
    const_cast<std::vector<double> &>(v_min_prediction).assign(size, value);
  }
  void set_v_max_prediction(std::size_t size, double value) {
    const_cast<std::vector<double> &>(v_max_prediction).assign(size, value);
  }
  void set_driving_cycle(std::vector<double> value) {
    const_cast<std::vector<double> &>(driving_cycle) = value;
  }
  void set_driving_cycle(std::vector<double> value, double Ts) {
    const_cast<std::vector<double> &>(driving_cycle) = value;
    const_cast<double &>(Ts) = Ts;
  }
  void set_vp_prediction(std::vector<double> value) {
    const_cast<std::vector<double> &>(vp_prediction) = value;
  }
  void set_alpha_prediction(std::vector<double> value) {
    const_cast<std::vector<double> &>(alpha_prediction) = value;
  }
  void set_v_min_prediction(std::vector<double> value) {
    const_cast<std::vector<double> &>(v_min_prediction) = value;
  }
  void set_v_max_prediction(std::vector<double> value) {
    const_cast<std::vector<double> &>(v_max_prediction) = value;
  }
  void set_alpha(double value) { const_cast<double &>(alpha) = value; }
  void set_v_min(double value) { const_cast<double &>(v_min) = value; }
  void set_v_max(double value) { const_cast<double &>(v_max) = value; }
  void set_f_t_prev(double value) { const_cast<double &>(f_t_prev) = value; }
  void set_f_b_prev(double value) { const_cast<double &>(f_b_prev) = value; }
  void set_v_p_prev(double value) { const_cast<double &>(v_p_prev) = value; }
  void set_v_h_prev(double value) { const_cast<double &>(v_h_prev) = value; }
  void set_v_h_0(double value) { const_cast<double &>(v_h_0) = value; }
  void set_d_h_0(double value) { const_cast<double &>(d_h_0) = value; }
};

} // namespace cg2o::mpc

#endif // MPC_PARAMETERS_H
