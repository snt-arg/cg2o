#ifndef G2O_BASE_FIXED_SIZED_EDGE_EQ
#define G2O_BASE_FIXED_SIZED_EDGE_EQ

#include "g2o/core/base_fixed_sized_edge.h"
#include "vertex_lagrange_multiplier.h"
#include <Eigen/Core>
#include <utility>
#include <vector>

namespace cg2o {

template <int D, typename E, typename... VertexTypes>
class G2O_CORE_API BaseFixedSizedEdgeEq
    : public g2o::BaseFixedSizedEdge<2 * D, E, VertexTypes...,
                                     VertexLagrangeMultiplier<D>> {
public:
  using VoidEdgeFuncType =
      std::function<void(BaseFixedSizedEdgeEq<D, E, VertexTypes...> &)>;
  using InformationType =
      typename g2o::BaseFixedSizedEdge<D, E, VertexTypes...>::InformationType;

  BaseFixedSizedEdgeEq();

  // User-defined equality computation
  virtual void computeEq() = 0;

  virtual double chi2() const override;

  void computeError() override;
  void constructQuadraticForm() override;

  void setVertexLagrangeMultiplierId(int id);
  int getVertexLagrangeMultiplierId() const;

  // Setter for the function pointers
  void setConstructQuadraticFormImpl(const VoidEdgeFuncType &func);

  /**
   * Dual (multiplier) regularization of this equality constraint.
   * The error is e = [h(x); nu] with information Omega = [[0, I], [I, 0]].
   * With sigma > 0 the Hessian block is built with
   *   Omega_H = [[0, I], [I, -sigma I]]
   * and the gradient with the exact Omega, so this constraint contributes
   *   [[0, J^T], [J, -sigma I]] to the Hessian (rows/cols x, nu) and
   *   -[J^T nu; h] to the right-hand side.
   * The sign must be negative (proximal method of multipliers); +sigma gives
   * a wrong step. sigma < 0 (default): use the optimizer's
   * setEqualityRegularization() value; sigma = 0: exact KKT for this edge.
   */
  void setDualRegularization(double sigma) { _dualRegularization = sigma; }
  double dualRegularization() const { return _dualRegularization; }
  //! assembles this edge's Hessian blocks and gradient with the given sigma
  void constructQuadraticFormRegularized(double sigma);
  /**
   * Closed-form assembly for the standard information matrix Omega =
   * [[0, I], [I, 0]] (ISPD, BIPM). With J_i the first D rows of the Jacobian
   * of vertex i (dh/dx_i) and nu the multiplier, the contributions are
   *   H(x_i, x_j) = 0,  H(x_i, nu) = J_i^T,  H(nu, nu) = -sigma I,
   *   b(x_i) = -J_i^T nu,  b(nu) = -h,
   * so the dense J^T Omega J products of the generic path are not needed.
   * Used automatically by constructQuadraticFormRegularized().
   */
  void constructQuadraticFormFast(double sigma);
  //! false: always use the generic g2o assembly (default true; for A/B testing)
  void setClosedFormAssembly(bool on) { _closedForm = on; }
  bool closedFormAssembly() const { return _closedForm; }

  // Setter
  void
  setLagrangeMultiplier(const Eigen::Matrix<double, D, 1> &lagrangeMuliplier);
  void setRho(const Eigen::Matrix<double, D, 1> &rho);
  void setConstraintViolationPrev(
      Eigen::Matrix<double, D, 1> &&input,
      size_t start = 0); // in AL method we use the information matrix to save
                         // the constraint violation and rho_bar

  // Getter functions
  Eigen::Matrix<double, D, 1> lagrangeMultiplier() const;
  Eigen::Matrix<double, D, 1> rho() const;

  bool read(std::istream &is) override;
  bool write(std::ostream &os) const override;
  void initializeInformationMatrix();

protected:
  void setInformation(const InformationType &information);
  Eigen::Matrix<double, D, 1> _eq = Eigen::Matrix<double, D, 1>::Zero();
  int _VertexLagrangeMultiplierId = -1;
  void initializeJacobians();

  VoidEdgeFuncType _constructQuadraticFormImpl = nullptr;
  double _dualRegularization = -1.0; // < 0: use the optimizer's default
  bool _standardOmega = false; // true while _information == [[0,I],[I,0]]
  bool _closedForm = true;

private:
  template <std::size_t I>
  void assemblePrimal(const Eigen::Matrix<double, D, 1> &nu);

protected:

  // In case Rho is used
  Eigen::Matrix<double, D, 1> _rho;
};

} // namespace cg2o

namespace cg2o {

template <int D, typename E, typename... VertexTypes>
BaseFixedSizedEdgeEq<D, E, VertexTypes...>::BaseFixedSizedEdgeEq() {
  this->_information.setZero();
}

template <int D, typename E, typename... VertexTypes>
double BaseFixedSizedEdgeEq<D, E, VertexTypes...>::chi2() const {
  return 0.0;
}
template <int D, typename E, typename... VertexTypes>
void BaseFixedSizedEdgeEq<D, E, VertexTypes...>::initializeJacobians() {
  // Helper lambda to set each Jacobian to zero
  auto setZeroHelper = [](auto &jacobian) { jacobian.setZero(); };

  // Apply the lambda to each Jacobian in the tuple
  std::apply(
      [&setZeroHelper](auto &...jacobians) {
        (setZeroHelper(jacobians),
         ...); // Expand and apply to all tuple elements
      },
      this->_jacobianOplus);

  // Set the Jacobian for the Lagrangian vertex
  if constexpr (sizeof...(VertexTypes) > 0) {
    auto &J_vLagrangian =
        std::get<sizeof...(VertexTypes)>(this->_jacobianOplus);
    J_vLagrangian.setZero();
    J_vLagrangian.template block<D, D>(D, 0).setIdentity();
  }
}

template <int D, typename E, typename... VertexTypes>
void BaseFixedSizedEdgeEq<D, E, VertexTypes...>::computeError() {
  if (this->_vertices.empty() ||
      this->_vertices.size() < sizeof...(VertexTypes) + 1) {
    throw std::runtime_error("Lagrangian vertex missing or invalid index");
  }

  const auto *vLagrangian =
      static_cast<const VertexLagrangeMultiplier<D> *>(this->_vertices.back());
  Eigen::Matrix<double, D, 1> Lagrangian = vLagrangian->estimate();

  // Compute equality using the user-defined method
  this->computeEq();

  // Update the error vector
  this->_error.segment(0, D) = _eq;
  this->_error.segment(D, D) = Lagrangian.segment(0, D);

#if CG2O_DEBUG_LINEAR_SOLVER
  if (!this->error().allFinite()) {
    std::cout << "========================================== " << std::endl;
    std::cout << "[WARNING] ****************** Invalid value detected in the "
                 "eq edge  "
                 "****************** [WARNING]"
              << std::endl;
    std::cout << "Edge Type: " << typeid(*this).name() << std::endl;
  }
  if (false) {
    std::cout << "Edge ID: " << this->id() << std::endl;
    std::cout << "Error vector Eq: " << this->error().transpose() << std::endl;
    std::cout << "Lagrange Multiplier Eq: "
              << this->lagrangeMultiplier().transpose() << std::endl;
    std::cout << "Rho Eq: " << this->rho().transpose() << std::endl;
    std::cout << "========================================== " << std::endl;
  }
#endif
}

template <int D, typename E, typename... VertexTypes>
void BaseFixedSizedEdgeEq<D, E, VertexTypes...>::initializeInformationMatrix() {
  this->_information.setZero();
  this->_information.template block<D, D>(0, D).setIdentity();
  this->_information.template block<D, D>(D, 0).setIdentity();
  _standardOmega = true;
}

template <int D, typename E, typename... VertexTypes>
int BaseFixedSizedEdgeEq<D, E, VertexTypes...>::getVertexLagrangeMultiplierId()
    const {
  return _VertexLagrangeMultiplierId;
}

template <int D, typename E, typename... VertexTypes>
bool BaseFixedSizedEdgeEq<D, E, VertexTypes...>::write(std::ostream &os) const {
  for (int i = 0; i < D; ++i) {
    os << this->_error[i] << " ";
  }
  return os.good();
}

template <int D, typename E, typename... VertexTypes>
bool BaseFixedSizedEdgeEq<D, E, VertexTypes...>::read(std::istream &is) {
  for (int i = 0; i < D; ++i) {
    is >> this->_error[i];
  }
  return is.good();
}

template <int D, typename E, typename... VertexTypes>
void BaseFixedSizedEdgeEq<D, E, VertexTypes...>::setConstructQuadraticFormImpl(
    const VoidEdgeFuncType &func) {
  _constructQuadraticFormImpl = func;
}

// Override the constructQuadraticForm method
template <int D, typename E, typename... VertexTypes>
void BaseFixedSizedEdgeEq<D, E, VertexTypes...>::constructQuadraticForm() {
  // Ensure that a function has been assigned before calling it
  if (_constructQuadraticFormImpl) {
    _constructQuadraticFormImpl(
        *this); // Call the assigned function (or lambda)
  } else if (_dualRegularization >= 0.0) {  // explicit sigma, 0 = exact KKT (closed form)
    constructQuadraticFormRegularized(_dualRegularization);
  } else {
    // If no function is assigned, call the base class implementation (if
    // needed)
    g2o::BaseFixedSizedEdge<
        2 * D, E, VertexTypes...,
        VertexLagrangeMultiplier<D>>::constructQuadraticForm();
  }
}

template <int D, typename E, typename... VertexTypes>
void BaseFixedSizedEdgeEq<D, E, VertexTypes...>::
    constructQuadraticFormRegularized(double sigma) {
  if (_standardOmega && _closedForm) {
    constructQuadraticFormFast(sigma);
    return;
  }
  constexpr std::size_t nv = sizeof...(VertexTypes) + 1; // + multiplier vertex
  // right-hand side with the exact information: -Omega e = -[nu; h]
  const auto weightedError = (-(this->_information * this->_error)).eval();
  if (sigma > 0.0) {
    auto omegaH = this->_information; // [[0, I], [I, 0]]
    omegaH.diagonal().template tail<D>().array() -= sigma;
    this->constructQuadraticFormNs(omegaH, weightedError,
                                   std::make_index_sequence<nv>());
  } else {
    this->constructQuadraticFormNs(this->_information, weightedError,
                                   std::make_index_sequence<nv>());
  }
}

template <int D, typename E, typename... VertexTypes>
void BaseFixedSizedEdgeEq<D, E, VertexTypes...>::constructQuadraticFormFast(
    double sigma) {
  constexpr std::size_t NV = sizeof...(VertexTypes); // multiplier vertex index
  const Eigen::Matrix<double, D, 1> h = this->_error.template head<D>();
  const Eigen::Matrix<double, D, 1> nu = this->_error.template tail<D>();

  auto *vnu = this->template vertexXn<NV>();
  if (!vnu->fixed()) {
    g2o::internal::QuadraticFormLock lck(*vnu);
    (void)lck;
    vnu->b().noalias() -= h;
    if (sigma > 0.0)
      vnu->A().diagonal().array() -= sigma; // H(nu, nu) = -sigma I
  }
  [&]<std::size_t... I>(std::index_sequence<I...>) {
    (this->template assemblePrimal<I>(nu), ...);
  }(std::make_index_sequence<NV>());
}

template <int D, typename E, typename... VertexTypes>
template <std::size_t I>
void BaseFixedSizedEdgeEq<D, E, VertexTypes...>::assemblePrimal(
    const Eigen::Matrix<double, D, 1> &nu) {
  constexpr std::size_t NV = sizeof...(VertexTypes);
  auto *v = this->template vertexXn<I>();
  if (v->fixed())
    return;
  const auto J = std::get<I>(this->_jacobianOplus).template topRows<D>();
  {
    g2o::internal::QuadraticFormLock lck(*v);
    (void)lck;
    v->b().noalias() -= J.transpose() * nu; // b(x_i) -= J_i^T nu
  }
  auto *vnu = this->template vertexXn<NV>();
  if (vnu->fixed())
    return;
  constexpr int K =
      g2o::internal::pair_to_index(static_cast<int>(I), static_cast<int>(NV));
  g2o::internal::QuadraticFormLock lck(*vnu);
  (void)lck;
  if (this->_hessianRowMajor[K])
    std::get<K>(this->_hessianTupleTransposed).noalias() += J;
  else
    std::get<K>(this->_hessianTuple).noalias() += J.transpose();
}

template <int D, typename E, typename... VertexTypes>
void BaseFixedSizedEdgeEq<D, E, VertexTypes...>::setInformation(
    const InformationType &information) {
  this->_information = information;
  _standardOmega = false;
}

template <int D, typename E, typename... VertexTypes>
void BaseFixedSizedEdgeEq<D, E, VertexTypes...>::setConstraintViolationPrev(
    Eigen::Matrix<double, D, 1> &&input, size_t start) {
  this->_information.diagonal().segment(start, D) = input;
  _standardOmega = false;
}

template <int D, typename E, typename... VertexTypes>
void BaseFixedSizedEdgeEq<D, E, VertexTypes...>::setVertexLagrangeMultiplierId(
    int id) {
  _VertexLagrangeMultiplierId = id;
}

// Setter functions
template <int D, typename E, typename... VertexTypes>
void BaseFixedSizedEdgeEq<D, E, VertexTypes...>::setLagrangeMultiplier(
    const Eigen::Matrix<double, D, 1> &lagrangeMultiplier) {
  // Ensure the last vertex is of type VertexLagrangeMultiplier<D>
  if (this->_vertices.empty() ||
      this->_vertices.size() < sizeof...(VertexTypes) + 1) {
    throw std::runtime_error("Lagrangian vertex missing or invalid index");
  }

  // Get the Lagrangian multiplier vertex
  auto *vLagrangian =
      static_cast<VertexLagrangeMultiplier<D> *>(this->_vertices.back());

  // Set the estimate of the Lagrangian multiplier vertex
  vLagrangian->setEstimate(lagrangeMultiplier);
}

template <int D, typename E, typename... VertexTypes>
void BaseFixedSizedEdgeEq<D, E, VertexTypes...>::setRho(
    const Eigen::Matrix<double, D, 1> &rho) {
  _rho = rho;
}

// Getter functions

template <int D, typename E, typename... VertexTypes>
Eigen::Matrix<double, D, 1>
BaseFixedSizedEdgeEq<D, E, VertexTypes...>::lagrangeMultiplier() const {
  // Ensure the vertex list is valid
  if (this->_vertices.empty() ||
      this->_vertices.size() < sizeof...(VertexTypes) + 1) {
    throw std::runtime_error("Lagrangian vertex missing or invalid index");
  }

  // Get the Lagrangian multiplier vertex
  const auto *vLagrangian =
      static_cast<const VertexLagrangeMultiplier<D> *>(this->_vertices.back());

  // Return the estimate of the Lagrange multiplier vertex
  return vLagrangian->estimate();
}

template <int D, typename E, typename... VertexTypes>
Eigen::Matrix<double, D, 1>
BaseFixedSizedEdgeEq<D, E, VertexTypes...>::rho() const {
  return _rho;
}

} // namespace cg2o

#endif // G2O_BASE_FIXED_SIZED_EDGE_EQ