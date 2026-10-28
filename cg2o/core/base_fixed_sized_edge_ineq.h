#ifndef G2O_BASE_FIXED_SIZED_EDGE_INEQ_H
#define G2O_BASE_FIXED_SIZED_EDGE_INEQ_H

#include "g2o/core/base_fixed_sized_edge.h"
#include <functional> // Required for std::function
#include <property.h>

namespace cg2o {

// Declare the template class
template <int D, typename E, typename... VertexTypes>
class G2O_CORE_API BaseFixedSizedEdgeIneq
    : public g2o::BaseFixedSizedEdge<D, E, VertexTypes...> {
public:
  // Type alias for the function pointer using std::function
  using VoidEdgeFuncType =
      std::function<void(BaseFixedSizedEdgeIneq<D, E, VertexTypes...> &)>;
  using InformationType =
      typename g2o::BaseFixedSizedEdge<D, E, VertexTypes...>::InformationType;

  BaseFixedSizedEdgeIneq();

  // Virtual method to compute inequality; must be implemented by derived
  // classes
  virtual void computeIneq() = 0;

  virtual double chi2() const override;

  // Override methods
  void computeError() override;
  void constructQuadraticForm() override;
  bool read(std::istream &is) override;
  bool write(std::ostream &os) const override;

  // Setter for the function pointers
  void setConstructQuadraticFormImpl(const VoidEdgeFuncType &func);

  // Setter
  void
  setLagrangeMultiplier(const Eigen::Matrix<double, D, 1> &lagrangeMultiplier);
  void setSlackVariable(const Eigen::Matrix<double, D, 1> &slackVariable);
  void setRho(const Eigen::Matrix<double, D, 1> &rho);

  auto &jacobianOplus() const { return this->_jacobianOplus; }
  void computeJacobianUpdateProduct(const double *update,
                                    Eigen::Matrix<double, D, 1> &result) const;

  // Getter functions
  Eigen::Matrix<double, D, 1> &lagrangeMultiplier();
  Eigen::Matrix<double, D, 1> &multiplierUpdate();
  Eigen::Matrix<double, D, 1> &slackVariable();
  Eigen::Matrix<double, D, 1> &slackUpdate();
  Eigen::Matrix<double, D, 1> &rho();

  void setConstraintViolationPrev(
      Eigen::Matrix<double, D, 1> &&input,
      size_t start = 0); // in AL method we use the information matrix to save
                         // the constraint violation

  // ---- warm start (used by SparseOptimizerISPD) ----
  // The optimizer keeps one stored copy of (lambda, slack, g(x)) per edge:
  // the iterate it selected as warm-start point (final or middle iterate of
  // the previous solve). `stamp` identifies the solve that stored it, so an
  // edge that was not active in that solve is recognised and initialised
  // cold. storeAuxiliary() expects computeError() to be current.
  void storeAuxiliary(long stamp) {
    _lagrangeMultiplierStored = _lagrangeMultiplier;
    _slackVariableStored = _slackVariable;
    _ineqStored = this->_error;
    _auxiliaryStamp = stamp;
  }
  //! copies the stored point back; false (nothing changed) if it is not
  //! from solve `stamp`
  bool restoreAuxiliary(long stamp) {
    if (_auxiliaryStamp != stamp)
      return false;
    _lagrangeMultiplier = _lagrangeMultiplierStored;
    _slackVariable = _slackVariableStored;
    _ineqWarm = _ineqStored;
    return true;
  }
  long auxiliaryStamp() const { return _auxiliaryStamp; }
  //! true while (lambda, slack) hold warm values for the coming solve
  bool auxiliaryWarm() const { return _auxiliaryWarm; }
  void setAuxiliaryWarm(bool warm) { _auxiliaryWarm = warm; }
  //! g(x) of the restored warm point (shifted with lambda and slack); the
  //! optimizer compares it with g at the start of the new solve
  Eigen::Matrix<double, D, 1> &warmIneq() { return _ineqWarm; }
  // Separate rollback state: selecting a middle point must not overwrite it.
  void storeStartingAuxiliary(long stamp) {
    _startMultiplier = _lagrangeMultiplier;
    _startSlack = _slackVariable;
    _startStamp = stamp;
  }
  bool restoreStartingAuxiliary(long stamp) {
    if (_startStamp != stamp)
      return false;
    _lagrangeMultiplier = _startMultiplier;
    _slackVariable = _startSlack;
    return true;
  }

protected:
  void setInformation(const InformationType &information);
  // Member to store the function pointer (or lambda)
  typename g2o::BaseFixedSizedEdge<D, E, VertexTypes...>::ErrorVector
      _ineq; // Vector to store the inequality values

  VoidEdgeFuncType _constructQuadraticFormImpl = nullptr;

  Eigen::Matrix<double, D, 1> _lagrangeMultiplier;
  Eigen::Matrix<double, D, 1> _multiplierUpdate;
  double _lagrangeMultiplierMaxUpdateSize;

  // In case the slack variable is used
  Eigen::Matrix<double, D, 1> _slackVariable;
  Eigen::Matrix<double, D, 1> _slackUpdate;
  double _slackVariableMaxUpdateSize;

  // In case Rho is used
  Eigen::Matrix<double, D, 1> _rho;

  // warm-start storage (see storeAuxiliary)
  Eigen::Matrix<double, D, 1> _lagrangeMultiplierStored;
  Eigen::Matrix<double, D, 1> _slackVariableStored;
  Eigen::Matrix<double, D, 1> _ineqStored;
  Eigen::Matrix<double, D, 1> _ineqWarm;
  long _auxiliaryStamp = -1;
  bool _auxiliaryWarm = false;
  Eigen::Matrix<double, D, 1> _startMultiplier, _startSlack;
  long _startStamp = -1;
};

// Forward declare the customQuadraticForm function template
template <int D, typename E, typename... VertexTypes>
void customQuadraticForm(BaseFixedSizedEdgeIneq<D, E, VertexTypes...> &edge);

} // namespace cg2o

namespace cg2o {

template <int D, typename E, typename... VertexTypes>
BaseFixedSizedEdgeIneq<D, E, VertexTypes...>::BaseFixedSizedEdgeIneq()
    : g2o::BaseFixedSizedEdge<D, E, VertexTypes...>(),
      _constructQuadraticFormImpl(nullptr) {
  this->_information.setZero();
  // Constructor body (if needed)
}

// Override the constructQuadraticForm method
template <int D, typename E, typename... VertexTypes>
void BaseFixedSizedEdgeIneq<D, E, VertexTypes...>::constructQuadraticForm() {
  // Ensure that a function has been assigned before calling it
  if (_constructQuadraticFormImpl) {
    _constructQuadraticFormImpl(
        *this); // Call the assigned function (or lambda)
  } else {
    // If no function is assigned, call the base class implementation (if
    // needed)
    g2o::BaseFixedSizedEdge<D, E, VertexTypes...>::constructQuadraticForm();
  }
}

template <int D, typename E, typename... VertexTypes>
void BaseFixedSizedEdgeIneq<D, E, VertexTypes...>::computeError() {
  this->computeIneq();
  this->_error = _ineq;

#if CG2O_DEBUG_LINEAR_SOLVER
  if (!this->error().allFinite()) {
    std::cout << "========================================== " << std::endl;
    std::cout << "[WARNING] ****************** Invalid value detected in the "
                 "ineq edge  "
                 "****************** [WARNING]"
              << std::endl;
    std::cout << "Edge Type: " << typeid(*this).name() << std::endl;
  }
  if (false) {
    std::cout << "Edge ID: " << this->id() << std::endl;
    std::cout << "Edge Type: " << typeid(*this).name() << std::endl;
    std::cout << "Error vector Ineq: " << this->error().transpose()
              << std::endl;
    std::cout << "Lagrange Multiplier Ineq: "
              << this->lagrangeMultiplier().transpose() << std::endl;
    std::cout << "Slack Variable: " << this->slackVariable().transpose()
              << std::endl;
    std::cout << "Rho: " << this->rho().transpose() << std::endl;
    std::cout << "========================================== " << std::endl;
  }
#endif
}

template <int D, typename E, typename... VertexTypes>
double BaseFixedSizedEdgeIneq<D, E, VertexTypes...>::chi2() const {
  return 0.0;
}

template <int D, typename E, typename... VertexTypes>
bool BaseFixedSizedEdgeIneq<D, E, VertexTypes...>::write(
    std::ostream &os) const {
  for (int i = 0; i < D; ++i) {
    os << this->_error[i] << " ";
  }
  return os.good();
}

template <int D, typename E, typename... VertexTypes>
bool BaseFixedSizedEdgeIneq<D, E, VertexTypes...>::read(std::istream &is) {
  for (int i = 0; i < D; ++i) {
    is >> this->_error[i];
  }
  return is.good();
}

template <int D, typename E, typename... VertexTypes>
void BaseFixedSizedEdgeIneq<D, E, VertexTypes...>::computeJacobianUpdateProduct(
    const double *update, Eigen::Matrix<double, D, 1> &result) const {
  auto processVertex = [&](size_t i, const auto &J_v_i) {
    const auto *vertex =
        static_cast<const g2o::OptimizableGraph::Vertex *>(this->vertex(i));

    if (!vertex || vertex->fixed()) {
      return;
    }

    const int offset = vertex->colInHessian();
    const int dim = vertex->dimension();

    Eigen::Map<const Eigen::VectorXd> updateVec(update + offset, dim);
    result += J_v_i * updateVec;
  };

  // Process each vertex with its corresponding Jacobian
  size_t i = 0;
  std::apply([&](const auto &...Js) { (processVertex(i++, Js), ...); },
             this->_jacobianOplus);
};

// Setter functions pointers
template <int D, typename E, typename... VertexTypes>
void BaseFixedSizedEdgeIneq<D, E, VertexTypes...>::setInformation(
    const InformationType &information) {
  this->_information = information;
}

template <int D, typename E, typename... VertexTypes>
void BaseFixedSizedEdgeIneq<D, E, VertexTypes...>::setConstraintViolationPrev(
    Eigen::Matrix<double, D, 1> &&input, size_t start) {
  this->_information.diagonal().segment(start, D) = input;
}

template <int D, typename E, typename... VertexTypes>
void BaseFixedSizedEdgeIneq<D, E, VertexTypes...>::
    setConstructQuadraticFormImpl(const VoidEdgeFuncType &func) {
  _constructQuadraticFormImpl = func;
}

// Setter functions
template <int D, typename E, typename... VertexTypes>
void BaseFixedSizedEdgeIneq<D, E, VertexTypes...>::setLagrangeMultiplier(
    const Eigen::Matrix<double, D, 1> &lagrangeMultiplier) {
  _lagrangeMultiplier = lagrangeMultiplier;
}

template <int D, typename E, typename... VertexTypes>
void BaseFixedSizedEdgeIneq<D, E, VertexTypes...>::setSlackVariable(
    const Eigen::Matrix<double, D, 1> &slackVariable) {
  _slackVariable = slackVariable;
}

template <int D, typename E, typename... VertexTypes>
void BaseFixedSizedEdgeIneq<D, E, VertexTypes...>::setRho(
    const Eigen::Matrix<double, D, 1> &rho) {
  _rho = rho;
}

// Getter functions

template <int D, typename E, typename... VertexTypes>
Eigen::Matrix<double, D, 1> &
BaseFixedSizedEdgeIneq<D, E, VertexTypes...>::lagrangeMultiplier() {
  return _lagrangeMultiplier;
}

template <int D, typename E, typename... VertexTypes>
Eigen::Matrix<double, D, 1> &
BaseFixedSizedEdgeIneq<D, E, VertexTypes...>::multiplierUpdate() {
  return _multiplierUpdate;
}

template <int D, typename E, typename... VertexTypes>
Eigen::Matrix<double, D, 1> &
BaseFixedSizedEdgeIneq<D, E, VertexTypes...>::slackVariable() {
  return _slackVariable;
}

template <int D, typename E, typename... VertexTypes>
Eigen::Matrix<double, D, 1> &
BaseFixedSizedEdgeIneq<D, E, VertexTypes...>::slackUpdate() {
  return _slackUpdate;
}

template <int D, typename E, typename... VertexTypes>
Eigen::Matrix<double, D, 1> &
BaseFixedSizedEdgeIneq<D, E, VertexTypes...>::rho() {
  return _rho;
}

} // namespace cg2o
#endif // G2O_BASE_FIXED_SIZED_EDGE_INEQ_H
