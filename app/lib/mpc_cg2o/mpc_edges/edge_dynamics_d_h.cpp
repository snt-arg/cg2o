#include "edge_dynamics_d_h.h"
#include <cmath>

namespace cg2o::mpc {

EdgeDynamics_d_h::EdgeDynamics_d_h(
    int k, std::shared_ptr<MPCParameters> param = nullptr)
    : _k(k), _param(param) {}

void EdgeDynamics_d_h::computeEq() {
  const Vertex_d_h *v_d_k = static_cast<const Vertex_d_h *>(this->_vertices[0]);
  const Vertex_d_h *v_d_kp1 =
      static_cast<const Vertex_d_h *>(this->_vertices[1]);
  const Vertex_v_h *v_v_k = static_cast<const Vertex_v_h *>(this->_vertices[2]);
  const Vertex_v_h *v_v_kp1 =
      static_cast<const Vertex_v_h *>(this->_vertices[3]);

  double v_p_k = _param->get_vp_prediction(
      _k); // the prediction of preceding velocity at time k
  double v_p_kp1 = _param->get_vp_prediction(_k + 1);
  double d_h_k = v_d_k->estimate();
  double d_h_kp1 = v_d_kp1->estimate();
  double v_h_k = v_v_k->estimate();
  double v_h_kp1 = v_v_kp1->estimate();
  
  if (_param->useSpaceDomain()) { // Space domain
    double delta_s = _param->get_delta_s();
    double delta_t = 2 * delta_s / (v_h_k + v_h_kp1);
    _eq[0] = d_h_kp1 - (d_h_k + delta_t / 2 * (v_p_k + v_p_kp1) - delta_s);
    return;
  }

  if (!_param->useSpaceDomain()) { // Time domain
    double delta_t = _param->get_delta_t();
    _eq[0] =
        d_h_kp1 - (d_h_k + delta_t / 2 * (v_p_k + v_p_kp1 - (v_h_k + v_h_kp1)));
  }
}
#ifndef MPC_USE_NUMERICAL_JACOBIAN
void EdgeDynamics_d_h::linearizeOplus() {

  initializeJacobians(); // handles bottomRows(D);
   
  const int D = _eq.size();
  auto &&J_v0 = std::get<0>(_jacobianOplus).topRows(D);       // d_h_k
  auto &&J_v1 = std::get<1>(_jacobianOplus).topRows(D);       // d_h_kp1
  auto &&J_v2 = std::get<2>(this->_jacobianOplus).topRows(D); // v_h_k
  auto &&J_v3 = std::get<3>(this->_jacobianOplus).topRows(D); // v_h_kp1
  J_v0 << -1.0;                                               // J(e_0,v_0[0])
  J_v1 << 1.0;                                                // J(e_0,v_0[1])

  double v_h_k =
      static_cast<const Vertex_v_h *>(this->_vertices[2])->estimate();

  if (_param->useSpaceDomain()) { // Space domain
    double v_h_kp1 =
        static_cast<const Vertex_v_h *>(this->_vertices[3])->estimate();
    double v_p_k = _param->get_vp_prediction(_k);
    double v_p_kp1 = _param->get_vp_prediction(_k + 1);
    double delta_s = _param->get_delta_s();
    double temp = -delta_s * (v_p_k + v_p_kp1);
    J_v2 << -temp / pow(v_h_k + v_h_kp1, 2); // J(e_0,v_0[2])
    J_v3 << -temp / pow(v_h_k + v_h_kp1, 2); // J(e_0,v_0[3])
    // _eq[0] =  -  delta_s  * (v_p_k + v_p_kp1) / (v_h_k + v_h_kp1);
    return;
  }

  if (!_param->useSpaceDomain()) { // Time domain
    double delta_t = _param->get_delta_t();
    // d_h_kp1 - (d_h_k + delta_t / 2 * (v_p_k + v_p_kp1 - (v_h_k + v_h_kp1)));
    J_v2 << delta_t / 2; // J(e_0,v_0[2])
    J_v3 << delta_t / 2; // J(e_0,v_0[3])
    return;
  }
}
#endif

bool EdgeDynamics_d_h::write(std::ostream &os) const {
  os << "EdgeDynamics_d_h, ";
  for (int i = 0; i < 1; i++) {
    os << _eq[i] << " ";
  }

  return os.good();
}

bool EdgeDynamics_d_h::read(std::istream &is) {
  for (int i = 0; i < 1; i++) {
    is >> _eq[i];
  }

  return is.good();
}

} // namespace cg2o::mpc
