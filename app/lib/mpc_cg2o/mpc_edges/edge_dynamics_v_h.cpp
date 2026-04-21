#include "edge_dynamics_v_h.h"
#include <cmath>

namespace cg2o::mpc {

EdgeDynamics_v_h::EdgeDynamics_v_h(
    int k, std::shared_ptr<MPCParameters> param = nullptr)
    : _k(k), _param(param) {}

void EdgeDynamics_v_h::computeEq() {
  const Vertex_v_h *v_v_k = static_cast<const Vertex_v_h *>(this->_vertices[0]);
  const Vertex_v_h *v_v_kp1 =
      static_cast<const Vertex_v_h *>(this->_vertices[1]);
  const Vertex_f_t *v_f_t_k =
      static_cast<const Vertex_f_t *>(this->_vertices[2]);
  const Vertex_f_b *v_f_b_k =
      static_cast<const Vertex_f_b *>(this->_vertices[3]);

  double v_h_k = v_v_k->estimate();
  double v_h_kp1 = v_v_kp1->estimate();
  double f_t_k = v_f_t_k->estimate();
  double f_b_k = v_f_b_k->estimate();

  double m_eq = _param->get_m_eq();
  double g = _param->get_g();
  double alpha_k = _param->get_alpha_prediction(_k);
  double c_r = _param->get_c_r();
  double c_a = _param->get_c_a();
 
                                             // domain
    if (_param->useSpaceDomain()) { // Space domain dynamics
    double m_total = _param->get_m_total();
    double delta_s = _param->get_delta_s();
    double A = 1 - 2 * c_a * delta_s / m_eq;
    double B = 2 * delta_s / m_eq;
    double w_1 = -B * m_total * g *
                 (sin(alpha_k * M_PI / 180) + c_r * cos(alpha_k * M_PI / 180));

    if (_param->isLinearDynamics()) { // linear in v_h_kp1
      _eq[0] = v_h_kp1 - sqrt(pow(v_h_k, 2) * (A) + B * (f_t_k - f_b_k) + w_1);
    } else {
      const Vertex_v_h *v_v_kp1 =
          static_cast<const Vertex_v_h *>(this->_vertices[1]);
      double v_h_kp1 = v_v_kp1->estimate();
      _eq[0] =
          pow(v_h_kp1, 2) - (pow(v_h_k, 2) * (A) + B * (f_t_k - f_b_k) + w_1);
    }

    return;
  }

  if (!_param->useSpaceDomain()) { // Time domain dynamics
    double f_roll = _param->get_c_r() * _param->get_m_total() * g *
                    cos(alpha_k * M_PI / 180);
    double f_grav = _param->get_m_total() * g * sin(alpha_k * M_PI / 180);
    double p1 = _param->get_p1();
    double p2 = _param->get_p2();
    double delta_t = _param->get_delta_t();
    double w_1 = -delta_t / m_eq * (f_grav + f_roll + p1);

    if (_param->isLinearDynamics()) {
      _eq[0] = v_h_kp1 - ((1 - delta_t / m_eq * p2) * v_h_k +
                          delta_t / m_eq * (f_t_k - f_b_k) + w_1);
    } else {
      _eq[0] =
          v_h_kp1 -
          (v_h_k + (delta_t / m_eq) *
                       ((f_t_k - f_b_k) -
                        (f_grav + f_roll + _param->get_c_a() * v_h_k * v_h_k)));
    }
  }
}

#if MPC_USE_NUMERICAL_JACOBIANS == 0
void EdgeDynamics_v_h::linearizeOplus() {
  // The derivative of _eq[0]  with respect to the vertex[0] (v_h_k)  is ...
  // [ -(1 - delta_t / m_eq * p2)] if linear dynamics  or [- (1 + delta_t / m_eq
  // * 2 * _param->get_c_a() * v_h_k)] if non-linear dynamics The derivative of
  // _eq[0]  with respect to the vertex[1] (v_h_kp1)  is 1 The derivative of
  // _eq[0]  with respect to the vertex[2] (f_t_k)  is  delta_t / m_eq The
  // derivative of _eq[0]  with respect to the vertex[3] (f_b_k)  is -delta_t /
  // m_eq Set the Jacobian's value for the scalar relationship we need here
  // carefully to define the Jacobian matrix because the solver might add
  // Langrangian varibales so to be general we define it in the equatily edges
  // as following the dimension of this initialize the Jacobians to zero

  initializeJacobians();
  // The Jacobian for the Lagrangian vertex
  double delta_t = _param->get_delta_t();
  double m_eq = _param->get_m_eq();
  double p2 = _param->get_p2();
 

  const Vertex_v_h *v_v_k = static_cast<const Vertex_v_h *>(this->_vertices[0]);
  double v_h_k = v_v_k->estimate();

  const int D = _eq.size();
  auto &&J_v0 = std::get<0>(_jacobianOplus).topRows(D);       // v_h_k
  auto &&J_v1 = std::get<1>(_jacobianOplus).topRows(D);       // v_h_kp1
  auto &&J_v2 = std::get<2>(this->_jacobianOplus).topRows(D); // f_t_k
  auto &&J_v3 = std::get<3>(this->_jacobianOplus).topRows(D); // f_b_k

  if (_param->useSpaceDomain()) { // Space domain
    double c_a = _param->get_c_a();
    double g = _param->get_g();
    double alpha_k = _param->get_alpha_prediction(_k);
    double c_r = _param->get_c_r();
    double m_total = _param->get_m_total();
    double delta_s = _param->get_delta_s();
    double A = 1 - 2 * c_a * delta_s / m_eq;
    double B = 2 * delta_s / m_eq;
    double w_1 = -B * m_total * g *
                 (sin(alpha_k * M_PI / 180) + c_r * cos(alpha_k * M_PI / 180));

    if (_param->isLinearDynamics()) {
      const Vertex_f_t *v_f_t_k =
          static_cast<const Vertex_f_t *>(this->_vertices[2]);
      const Vertex_f_b *v_f_b_k =
          static_cast<const Vertex_f_b *>(this->_vertices[3]);
      double f_t_k = v_f_t_k->estimate();
      double f_b_k = v_f_b_k->estimate();
      double S = pow(v_h_k, 2) * A + B * (f_t_k - f_b_k) + w_1;
      S = std::max(S, 1e-12);
      double temp = std::sqrt(S);

      J_v0 << -(v_h_k * A) / temp;
      J_v1 << 1.0;
      J_v2 << -(0.5 * B) / temp;
      J_v3 << +(0.5 * B) / temp; // later

    } else {
      const Vertex_v_h *v_v_kp1 =
          static_cast<const Vertex_v_h *>(this->_vertices[1]);
      double v_h_kp1 = v_v_kp1->estimate();
      // pow(v_h_kp1, 2) - (pow(v_h_k, 2) * (A) + B * (f_t_k - f_b_k) + w_1);

      J_v0 << -2 * v_h_k * A;
      J_v1 << 2 * v_h_kp1;
      J_v2 << -B;
      J_v3 << B;
    }
    return;
  }

  if (!_param->useSpaceDomain()) { // Time domain
    if (_param->isLinearDynamics()) {
      J_v0 << -(1 - delta_t / m_eq * p2);
    } else {
      J_v0 << -(1 - delta_t / m_eq * 2 * _param->get_c_a() * v_h_k);
    }
    J_v1 << 1;
    J_v2 << -delta_t / m_eq;
    J_v3 << delta_t / m_eq;
    return;
  }
}
#endif

bool EdgeDynamics_v_h::write(std::ostream &os) const {
  os << "EdgeDynamics_v_h, ";
  for (int i = 0; i < 1; i++) {
    os << _eq[i] << " ";
  }

  return os.good();
}

bool EdgeDynamics_v_h::read(std::istream &is) {
  for (int i = 0; i < 1; i++) {
    is >> _eq[i];
  }

  return is.good();
}

} // namespace cg2o::mpc
