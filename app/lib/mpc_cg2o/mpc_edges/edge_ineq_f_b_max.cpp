#include "edge_ineq_f_b_max.h"
#include <cmath>

namespace cg2o::mpc {

EdgeIneq_f_b_max::EdgeIneq_f_b_max(
    int k, std::shared_ptr<MPCParameters> param = nullptr)
    : _k(k), _param(param) {
  if (_param->isScaleInequalities()) {
    _scaling_factor = {10e-5, 10e-5};
  } else {
    _scaling_factor = {1.0, 1.0};
  }
}

void EdgeIneq_f_b_max::computeIneq() {
  const Vertex_f_b *v_f_b_k = static_cast<const Vertex_f_b *>(_vertices[0]);
  const Vertex_v_h *v_v_h_k = static_cast<const Vertex_v_h *>(_vertices[1]);

  double f_b_k = v_f_b_k->estimate();
  double v_h_k = v_v_h_k->estimate();

  double brems_a1 = _param->get_Brems_a1();
  double brems_b1 = _param->get_Brems_b1();
  double b_max = _param->get_b_max();
  double regen_max = _param->get_Regen_max();
  if (_param->isLinearInequalities() || true) {
    _ineq[0] =
        f_b_k -
        (b_max + regen_max); // g(x) <= 0 means f_b(k) <= b_max + Regen_max
    _ineq[1] = f_b_k - (brems_a1 * v_h_k + b_max + brems_b1); // g(x) <= 0
  } else {
    _ineq[0] =
        pow(f_b_k, 3) - pow(b_max + regen_max,
                            3); // g(x) <= 0 means f_b(k) <= b_max + Regen_max
    _ineq[1] = pow(f_b_k, 3) - pow(brems_a1 * v_h_k + b_max + brems_b1, 3);

    _ineq[0] *= _scaling_factor[0];
    _ineq[1] *= _scaling_factor[1];
  }
}
#if MPC_USE_NUMERICAL_JACOBIANS == 0
void EdgeIneq_f_b_max::linearizeOplus() {
  // The derivative of _ineq[0]  with respect to the vertex[0] (f_b_k)  is 1
  // The derivative of _ineq[1]  with respect to the vertex[0] (f_b_k)  is 1
  // The derivative of _ineq[0]  with respect to the vertex[1] (v_h_k)  is 0
  // The derivative of _ineq[1]  with respect to the vertex[1] (v_h_k)  is
  // -brems_a1 Set the Jacobian's value for the scalar relationship
  const Vertex_f_b *v_f_b_k = static_cast<const Vertex_f_b *>(_vertices[0]);
  const Vertex_v_h *v_v_h_k = static_cast<const Vertex_v_h *>(_vertices[1]);

  double f_b_k = v_f_b_k->estimate();
  double v_h_k = v_v_h_k->estimate();

  double brems_a1 = _param->get_Brems_a1();
  double brems_b1 = _param->get_Brems_b1();
  double b_max = _param->get_b_max();

  auto &J_v0 = std::get<0>(this->_jacobianOplus); // f_b_k
  auto &J_v1 = std::get<1>(this->_jacobianOplus); // v_h_k
  if (_param->isLinearInequalities() || true) {
    J_v0 << 1.0,   // J(e_0,v_0[0])
        1.0;       // J(e_1,v_0[0])
    J_v1 << 0.0,   // J(e_0,v_1[0])
        -brems_a1; // J(e_1,v_1[0])
  } else {
    J_v0 << 3.0 * pow(f_b_k, 2) * _scaling_factor[0], // J(e_0,v_0[0])
        3.0 * pow(f_b_k, 2) * _scaling_factor[1];     // J(e_1,v_0[0])

    J_v1 << 0.0, // J(e_0,v_1[0])
        -3 * _param->get_Brems_a1() *
            pow(brems_a1 * v_h_k + b_max + brems_b1, 2) * _scaling_factor[1];
    ; // J(e_1,v_1[0])
  }
}
#endif

bool EdgeIneq_f_b_max::write(std::ostream &os) const {
  os << "EdgeIneq_f_b_max: ";
  for (int i = 0; i < 2; ++i) {
    os << "ineq[" << i << "] = " << _ineq[i] << " ";
  }
  return os.good();
}

bool EdgeIneq_f_b_max::read(std::istream &is) {
  for (int i = 0; i < 2; ++i) {
    is >> _ineq[i];
  }
  return is.good();
}

} // namespace cg2o::mpc
