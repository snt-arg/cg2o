#include "edge_ineq_f_t_max.h"
#include <cmath>

namespace cg2o::mpc {

EdgeIneq_f_t_max::EdgeIneq_f_t_max(
    int k, std::shared_ptr<MPCParameters> param = nullptr)
    : _k(k), _param(param) {
  if (_param->isScaleInequalities()) {
    _scaling_factor = {10e-5, 10e-5};
  } else {
    _scaling_factor = {1.0, 1.0};
  }
}

void EdgeIneq_f_t_max::computeIneq() {
  const Vertex_f_t *v_f_t_k = static_cast<const Vertex_f_t *>(_vertices[0]);
  const Vertex_v_h *v_v_h_k = static_cast<const Vertex_v_h *>(_vertices[1]);

  double f_t_k = v_f_t_k->estimate();
  double v_h_k = v_v_h_k->estimate();

  double a2 = _param->get_a2();
  double b2 = _param->get_b2();
  double b4 = _param->get_b4();

  if (_param->isLinearInequalities() || true) {
    _ineq[0] = f_t_k - b4; // g(x) <= 0 means f_t(k) <= b4
    _ineq[1] =
        f_t_k - (a2 * v_h_k + b2); // g(x) <= 0 means f_t(k) <= a2*v_h(k) + b2
  } else {
    _ineq[0] = pow(f_t_k, 3) - pow(b4, 3); // g(x) <= 0 means f_t(k) <= b4
    _ineq[1] =
        pow(f_t_k, 3) -
        pow(a2 * v_h_k + b2, 3); // g(x) <= 0 means f_t(k) <= a2*v_h(k) + b2 }
    _ineq[0] *= _scaling_factor[0];
    _ineq[1] *= _scaling_factor[1];
  }


}
#ifdef USE_EXACT_JACOBIANS
void EdgeIneq_f_t_max::linearizeOplus() {
  // The derivative of _ineq[0]  with respect to the vertex[0] (f_t_k)
  // The derivative of _ineq[1]  with respect to the vertex[0] (f_t_k)
  // The derivative of _ineq[0]  with respect to the vertex[1] (v_h_k)
  // The derivative of _ineq[1]  with respect to the vertex[1] (v_h_k)
  // Set the Jacobian's value for the scalar relationship
  const Vertex_f_t *v_f_t_k = static_cast<const Vertex_f_t *>(_vertices[0]);
  const Vertex_v_h *v_v_h_k = static_cast<const Vertex_v_h *>(_vertices[1]);

  double f_t_k = v_f_t_k->estimate();
  double v_h_k = v_v_h_k->estimate();

  double a2 = _param->get_a2();
  double b2 = _param->get_b2();

  auto &J_v0 = std::get<0>(this->_jacobianOplus); // f_t_k
  auto &J_v1 = std::get<1>(this->_jacobianOplus); // v_h_k
  if (_param->isLinearInequalities() || true) {
    J_v0 << 1.0, // J(e_0,v_0[0])
        1.0;     // J(e_1,v_0[0])
    J_v1 << 0.0, // J(e_0,v_1[0])
        -a2;     // J(e_1,v_1[0])
  } else {
    J_v0 << 3.0 * pow(f_t_k, 2) * _scaling_factor[0], // J(e_0,v_0[0])
        3.0 * pow(f_t_k, 2) * _scaling_factor[1];     // J(e_1,v_0[0])
    J_v1 << 0.0,                                      // J(e_0,v_1[0])
        -3 * _param->get_a2() * pow(a2 * v_h_k + b2, 2) *
            _scaling_factor[1]; // J(e_1,v_1[0])
  }
}
#endif

bool EdgeIneq_f_t_max::write(std::ostream &os) const {
  os << "EdgeIneq_f_t_max: ";
  for (int i = 0; i < 2; ++i) {
    os << "ineq[" << i << "] = " << _ineq[i] << " ";
  }
  return os.good();
}

bool EdgeIneq_f_t_max::read(std::istream &is) {
  for (int i = 0; i < 2; ++i) {
    is >> _ineq[i];
  }
  return is.good();
}

} // namespace cg2o::mpc
