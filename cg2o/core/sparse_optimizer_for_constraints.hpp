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

#include "g2o/core/optimization_algorithm.h"
#include "g2o/core/optimization_algorithm_with_hessian.h"
#include "g2o/core/solver.h"
#include "g2o/stuff/logger.h"
#include "sparse_optimizer_for_constraints.h"
#include <Eigen/Core>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <unordered_set>
#include <utility>
#include <vector>

#ifndef NDEBUG
namespace {
/**
 * tests whether there is a NaN in the array
 */
bool arrayHasNaN(const double *array, int size, int *nanIndex = 0) {
  for (int i = 0; i < size; ++i)
    if (std::isnan(array[i])) {
      if (nanIndex)
        *nanIndex = i;
      return true;
    }
  return false;
}
} // namespace
#endif

namespace cg2o {
using namespace std;

// Default constructor
template <typename Derived>
SparseOptimizerForConstraints<Derived>::SparseOptimizerForConstraints() {
  defineAlphaBacktracking(_alpha);
};

// Default destructor
template <typename Derived>
SparseOptimizerForConstraints<Derived>::~SparseOptimizerForConstraints() =
    default;

template <typename Derived>
void SparseOptimizerForConstraints<Derived>::invalidateActiveStructure() {
  clearIndexMapping();
  _ivMap.clear();
  _activeVertices.clear();
  _activeEdges.clear();
  _activeEdgesEq.clear();
  _activeEdgesIneq.clear();
  _fixedEdgesEq.clear();
  _fixedEdgesIneq.clear();
}

template <typename Derived>
bool SparseOptimizerForConstraints<Derived>::fixedConstraintsFeasible() {
  return verifyEqFeasibility(_fixedEdgesEq, -1.0) &&
         verifyIneqFeasibility(_fixedEdgesIneq, -1.0);
}

template <typename Derived>
bool SparseOptimizerForConstraints<Derived>::removeEdge(
    HyperGraph::Edge *edge) {
  if (!edges().count(edge))
    return false;
  auto *e = static_cast<OptimizableGraph::Edge *>(edge);
  OptimizableGraph::Vertex *multiplier = nullptr; // created by addEdgeEq()
  std::size_t multiplierSlot = 0;
  if (_edgeEqSet.count(e) && !e->vertices().empty()) {
    multiplierSlot = e->vertices().size() - 1;
    multiplier =
        static_cast<OptimizableGraph::Vertex *>(e->vertices()[multiplierSlot]);
    std::erase(_vEqLagrangeMultipliers, multiplier);
  }
  _edgeEqSet.erase(e);
  _edgeIneqSet.erase(e);
  updateMultipliersEqFunctionMap.erase(e);
  updateMultipliersIneqFunctionMap.erase(e);
  // update the active structure in place: optimize() can follow directly,
  // without initializeOptimization()
  removeFromActiveStructure(e);
  if (!g2o::SparseOptimizer::removeEdge(edge))
    return false;
  // the multiplier goes with its edge, unless another edge still uses it or
  // it is the vertex being removed (g2o removes its edges first)
  if (multiplier && multiplier->edges().empty() &&
      std::find(_removing_vertices.begin(), _removing_vertices.end(),
                multiplier) == _removing_vertices.end() &&
      static_cast<Derived *>(this)->removeVertex(multiplier)) {
#if !G2O_DELETE_IMPLICITLY_OWNED_OBJECTS
    // the edge stays with the caller (with implicit ownership g2o deleted
    // it and the multiplier)
    e->setVertex(multiplierSlot, nullptr);
    delete multiplier;
#endif
  }
  return true;
}

template <typename Derived>
void SparseOptimizerForConstraints<Derived>::removeFromActiveStructure(
    OptimizableGraph::Edge *edge) {
  std::erase(_fixedEdgesEq, edge);
  std::erase(_fixedEdgesIneq, edge);
  // an inactive edge (fixed, other level, not initialized yet) does not
  // decide whether any vertex is active
  if (std::erase(_activeEdges, edge) == 0)
    return;
  std::erase(_activeEdgesIneq, edge);
  std::erase(_activeEdgesEq, edge);

  // initializeOptimization() keeps a vertex active only while an active edge
  // uses it; one pass over the remaining active edges finds the vertices of
  // the removed edge that are still used (e.g. not an equality multiplier)
  const auto &vertices = edge->vertices();
  std::vector<char> stillUsed(vertices.size(), 0);
  for (auto *other : _activeEdges)
    for (auto *v : other->vertices())
      for (std::size_t i = 0; i < vertices.size(); ++i)
        if (vertices[i] == v)
          stillUsed[i] = 1;

  bool indexMapChanged = false;
  for (std::size_t i = 0; i < vertices.size(); ++i) {
    auto *v = static_cast<OptimizableGraph::Vertex *>(vertices[i]);
    // erase() returns 0 for a vertex repeated in the edge (already dropped)
    if (stillUsed[i] || !v || std::erase(_activeVertices, v) == 0)
      continue;
    const int idx = v->hessianIndex();
    if (idx >= 0 && idx < static_cast<int>(_ivMap.size()) && _ivMap[idx] == v) {
      _ivMap.erase(_ivMap.begin() + idx);
      indexMapChanged = true;
    }
    v->setHessianIndex(-1);
  }
  // same consecutive numbering as buildIndexMapping() on the remaining set;
  // optimize() rebuilds the linear system structure from it
  if (indexMapChanged)
    for (std::size_t k = 0; k < _ivMap.size(); ++k)
      _ivMap[k]->setHessianIndex(static_cast<int>(k));
}

template <typename Derived>
bool SparseOptimizerForConstraints<Derived>::removeVertex(HyperGraph::Vertex *v,
                                                          bool detach) {
  if (!v || vertex(v->id()) != v)
    return false;
  std::erase(_vEqLagrangeMultipliers,
             static_cast<OptimizableGraph::Vertex *>(v));
  invalidateActiveStructure();
  _removing_vertices.push_back(v);
  const bool removed = g2o::SparseOptimizer::removeVertex(v, detach);
  _removing_vertices.pop_back();
  return removed;
}

template <typename Derived>
void SparseOptimizerForConstraints<Derived>::clear() {
  _edgeEqSet.clear();
  _edgeIneqSet.clear();
  _vEqLagrangeMultipliers.clear();
  updateMultipliersEqFunctionMap.clear();
  updateMultipliersIneqFunctionMap.clear();
  invalidateActiveStructure();
  g2o::SparseOptimizer::clear();
}

template <typename Derived>
bool SparseOptimizerForConstraints<Derived>::attachConstraintEdge(
    OptimizableGraph::Edge *edge) {
  if (!edge || edges().count(edge))
    return false;
  for (auto *v : edge->vertices())
    if (!v || vertex(v->id()) != v)
      return false;
  if (g2o::OptimizableGraph::addEdge(edge))
    return true;
  // g2o may insert an edge before parameter/cache resolution fails. Detach
  // without deleting: a failed insertion leaves ownership with the caller.
  if (edges().erase(edge))
    for (auto *v : edge->vertices())
      if (v)
        v->edges().erase(edge);
  return false;
}

template <typename Derived>
template <int D, typename E, typename... VertexTypes>
VertexLagrangeMultiplier<D> *
SparseOptimizerForConstraints<Derived>::attachEqualityEdge(
    BaseFixedSizedEdgeEq<D, E, VertexTypes...> *edge) {
  if (!edge || edges().count(edge))
    return nullptr;
  const int id = edge->getVertexLagrangeMultiplierId();
  if (id < 0 || vertex(id))
    return nullptr;
  for (std::size_t i = 0; i + 1 < edge->vertices().size(); ++i) {
    auto *v = edge->vertices()[i];
    if (!v || vertex(v->id()) != v)
      return nullptr;
  }
  auto multiplier = std::make_unique<VertexLagrangeMultiplier<D>>();
  multiplier->setId(id);
  multiplier->setInitialValue(_lagrange_multiplier_initial_eq);
  if (!static_cast<Derived *>(this)->addVertex(multiplier.get()))
    return nullptr;
  auto *nu = multiplier.release();
  edge->setVertex(edge->vertices().size() - 1, nu);
  if (attachConstraintEdge(edge))
    return nu;
  edge->setVertex(edge->vertices().size() - 1, nullptr);
  static_cast<Derived *>(this)->removeVertex(nu);
#if !G2O_DELETE_IMPLICITLY_OWNED_OBJECTS
  delete nu;
#endif
  return nullptr;
}

// Setter solver paramters
template <typename Derived>
bool SparseOptimizerForConstraints<Derived>::initializeOptimization(
    HyperGraph::VertexSet &vset, int level) {
  if (edges().size() == 0) {
    G2O_WARN("Attempt to initialize an empty graph");
    return false;
  }
  preIteration(-1);
  bool workspaceAllocated = _jacobianWorkspace.allocate();
  (void)workspaceAllocated;
  assert(workspaceAllocated &&
         "Error while allocating memory for the Jacobians");
  clearIndexMapping();
  _activeVertices.clear();
  _activeVertices.reserve(vset.size());
  _activeEdges.clear();
  // the constraint lists are rebuilt below; without clearing them a second
  // initializeOptimization() call listed every constraint edge twice, and
  // ISPD then updated its multiplier and slack twice per iteration
  _activeEdgesIneq.clear();
  _activeEdgesEq.clear();
  _fixedEdgesEq.clear();
  _fixedEdgesIneq.clear();
  set<Edge *> auxEdgeSet; // temporary structure to avoid duplicates
  set<Edge *> fixedEdgeSet;
  for (HyperGraph::VertexSet::iterator it = vset.begin(); it != vset.end();
       ++it) {
    OptimizableGraph::Vertex *v = (OptimizableGraph::Vertex *)*it;
    const OptimizableGraph::EdgeSet &vEdges = v->edges();
    // count if there are edges in that level. If not remove from the pool
    int levelEdges = 0;
    for (OptimizableGraph::EdgeSet::const_iterator it = vEdges.begin();
         it != vEdges.end(); ++it) {
      OptimizableGraph::Edge *e =
          reinterpret_cast<OptimizableGraph::Edge *>(*it);
      if (level < 0 || e->level() == level) {
        bool allVerticesOK = true;
        for (vector<HyperGraph::Vertex *>::const_iterator vit =
                 e->vertices().begin();
             vit != e->vertices().end(); ++vit) {
          if (vset.find(*vit) == vset.end()) {
            allVerticesOK = false;
            break;
          }
        }
        if (allVerticesOK) {
          if (!e->allVerticesFixed()) {
            auxEdgeSet.insert(e);
            levelEdges++;
          } else
            fixedEdgeSet.insert(e);
        }
      }
    }
    if (levelEdges) {
      _activeVertices.push_back(v);

      // test for NANs in the current estimate if we are debugging
#ifndef NDEBUG
      int estimateDim = v->estimateDimension();
      if (estimateDim > 0) {
        g2o::VectorX estimateData(estimateDim);
        if (v->getEstimateData(estimateData.data()) == true) {
          int k;
          bool hasNan = arrayHasNaN(estimateData.data(), estimateDim, &k);
          if (hasNan)
            G2O_WARN("Vertex {} contains a nan entry at index {}", v->id(), k);
        }
      }
#endif
    }
  }

  _activeEdges.reserve(auxEdgeSet.size());
  for (auto *e : fixedEdgeSet) {
    if (_edgeEqSet.count(e))
      _fixedEdgesEq.push_back(e);
    if (_edgeIneqSet.count(e))
      _fixedEdgesIneq.push_back(e);
  }
  _activeEdgesIneq.reserve(_edgeIneqSet.size());
  for (set<Edge *>::iterator it = auxEdgeSet.begin(); it != auxEdgeSet.end();
       ++it) {
    _activeEdges.push_back(*it);
    //-----------------------------------------
    // add the active edge (it) to the active inequality edges if it is the
    // inequality set
    if (_edgeIneqSet.find(*it) != _edgeIneqSet.end())
      // Element found, add it to _activeEdgesIneq
      _activeEdgesIneq.push_back(
          *it); // Dereference the iterator to get the element
    // add the active edge (it) to the active equality edges if it is the
    // equality set
    if (_edgeEqSet.find(*it) != _edgeEqSet.end())
      // Element found, add it to _activeEdgesIneq
      _activeEdgesEq.push_back(
          *it); // Dereference the iterator to get the element
    //-----------------------------------------
  }
  sortVectorContainers();
  bool indexMappingStatus = buildIndexMapping(_activeVertices);
  postIteration(-1);
  return indexMappingStatus;
}

template <typename Derived>
bool SparseOptimizerForConstraints<Derived>::initializeOptimization(int level) {
  HyperGraph::VertexSet vset;
  for (VertexIDMap::iterator it = vertices().begin(); it != vertices().end();
       ++it)
    vset.insert(it->second);
  return initializeOptimization(vset, level);
}

template <typename Derived>
bool SparseOptimizerForConstraints<Derived>::updateInitialization(
    HyperGraph::VertexSet &vset, HyperGraph::EdgeSet &eset) {
  if (!_algorithm)
    return false;
  // O(1) membership test for _ivMap; a stale Hessian index fails the check
  auto inIndexMap = [this](OptimizableGraph::Vertex *v) {
    const int idx = v->hessianIndex();
    return idx >= 0 && idx < static_cast<int>(_ivMap.size()) &&
           _ivMap[idx] == v;
  };
  // Validate the whole request before mutating an existing index map.
  for (auto *raw : vset) {
    auto *v = static_cast<OptimizableGraph::Vertex *>(raw);
    if (!v || vertex(v->id()) != v || v->marginalized())
      return false;
  }
  for (auto *raw : eset) {
    auto *e = static_cast<OptimizableGraph::Edge *>(raw);
    if (!e || !edges().count(e) || e->numUndefinedVertices())
      return false;
    for (auto *rawVertex : e->vertices()) {
      auto *v = static_cast<OptimizableGraph::Vertex *>(rawVertex);
      if (!v->fixed() && !inIndexMap(v) && !vset.count(v))
        return false;
    }
  }
  // hash sets instead of linear searches through the containers
  std::unordered_set<OptimizableGraph::Vertex *> activeVertices(
      _activeVertices.begin(), _activeVertices.end());
  std::unordered_set<OptimizableGraph::Edge *> activeEdges(_activeEdges.begin(),
                                                           _activeEdges.end());
  std::unordered_set<OptimizableGraph::Edge *> fixedEdgesEq(
      _fixedEdgesEq.begin(), _fixedEdgesEq.end());
  std::unordered_set<OptimizableGraph::Edge *> fixedEdgesIneq(
      _fixedEdgesIneq.begin(), _fixedEdgesIneq.end());
  std::vector<HyperGraph::Vertex *> newVertices;
  HyperGraph::EdgeSet newEdges;
  for (auto *raw : vset) {
    auto *v = static_cast<OptimizableGraph::Vertex *>(raw);
    if (activeVertices.insert(v).second)
      _activeVertices.push_back(v);
    if (v->fixed())
      v->setHessianIndex(-1);
    else if (!inIndexMap(v)) {
      v->setHessianIndex(static_cast<int>(_ivMap.size()));
      _ivMap.push_back(v);
      newVertices.push_back(v);
    }
  }
  for (auto *raw : eset) {
    auto *e = static_cast<OptimizableGraph::Edge *>(raw);
    if (e->allVerticesFixed()) {
      if (_edgeEqSet.count(e) && fixedEdgesEq.insert(e).second)
        _fixedEdgesEq.push_back(e);
      if (_edgeIneqSet.count(e) && fixedEdgesIneq.insert(e).second)
        _fixedEdgesIneq.push_back(e);
      continue;
    }
    if (!activeEdges.insert(e).second)
      continue;
    _activeEdges.push_back(e);
    newEdges.insert(e);
    if (_edgeIneqSet.count(e))
      _activeEdgesIneq.push_back(e);
    if (_edgeEqSet.count(e))
      _activeEdgesEq.push_back(e);
  }
  sortVectorContainers();
  if (newVertices.empty() && newEdges.empty())
    return true;
  // addEdge() updates maximum workspace sizes, but does not allocate them.
  // New constraint dimensions/arity must be allocated before the next solve.
  if (!_jacobianWorkspace.allocate() ||
      !_algorithm->updateStructure(newVertices, newEdges)) {
    invalidateActiveStructure();
    return false;
  }
  return true;
}

template <typename Derived>
void SparseOptimizerForConstraints<Derived>::updateEdgeVertices(
    OptimizableGraph::Edge *edge, const double *update, double scalingFactor) {
  const auto &vertices = edge->vertices();
  for (std::size_t i = 0; i < vertices.size(); ++i) {
    auto *vertex = static_cast<OptimizableGraph::Vertex *>(vertices[i]);
    if (vertex->fixed() || std::find(vertices.begin(), vertices.begin() + i,
                                     vertex) != vertices.begin() + i)
      continue;
    const double *increment = update + vertex->colInHessian();
    const int dimension = vertex->dimension();
    if (scalingFactor != 1.0) {
      if (_scaledUpdate.size() < dimension)
        _scaledUpdate.resize(dimension);
      _scaledUpdate.head(dimension) =
          scalingFactor *
          Eigen::Map<const Eigen::VectorXd>(increment, dimension);
      increment = _scaledUpdate.data();
    }
    vertex->oplus(increment);
  }
}

template <typename Derived>
void SparseOptimizerForConstraints<Derived>::applyScaledUpdate(
    const double *update, double scale) {
  if (scale == 1.0) {
    g2o::SparseOptimizer::update(update);
    return;
  }
  const int size =
      static_cast<g2o::OptimizationAlgorithmWithHessian *>(_algorithm)
          ->solver()
          .vectorSize();
  if (_scaledUpdate.size() < size)
    _scaledUpdate.resize(size);
  _scaledUpdate.head(size) =
      scale * Eigen::Map<const Eigen::VectorXd>(update, size);
  g2o::SparseOptimizer::update(_scaledUpdate.data());
}

template <typename Derived>
bool SparseOptimizerForConstraints<Derived>::initializeOptimization(
    HyperGraph::EdgeSet &eset) {
  preIteration(-1);
  bool workspaceAllocated = _jacobianWorkspace.allocate();
  (void)workspaceAllocated;
  assert(workspaceAllocated &&
         "Error while allocating memory for the Jacobians");
  clearIndexMapping();
  _activeVertices.clear();
  _activeEdges.clear();
  _activeEdges.reserve(eset.size());
  //-----------------------------------------
  _activeEdgesIneq.clear();
  _activeEdgesEq.clear();
  _fixedEdgesEq.clear();
  _fixedEdgesIneq.clear();
  _activeEdgesIneq.reserve(_edgeIneqSet.size());
  //-----------------------------------------
  set<Vertex *> auxVertexSet; // temporary structure to avoid duplicates
  for (HyperGraph::EdgeSet::iterator it = eset.begin(); it != eset.end();
       ++it) {
    OptimizableGraph::Edge *e = (OptimizableGraph::Edge *)(*it);
    if (e->numUndefinedVertices())
      continue;
    if (e->allVerticesFixed()) {
      if (_edgeEqSet.count(e))
        _fixedEdgesEq.push_back(e);
      if (_edgeIneqSet.count(e))
        _fixedEdgesIneq.push_back(e);
      continue;
    }
    for (vector<HyperGraph::Vertex *>::const_iterator vit =
             e->vertices().begin();
         vit != e->vertices().end(); ++vit) {
      auxVertexSet.insert(static_cast<OptimizableGraph::Vertex *>(*vit));
    }
    _activeEdges.push_back(reinterpret_cast<OptimizableGraph::Edge *>(*it));
    //-----------------------------------------
    // add the inequality edges to the active edges if they are in the
    // inequality set
    if (_edgeIneqSet.find(static_cast<OptimizableGraph::Edge *>(*it)) !=
        _edgeIneqSet.end())
      // Element found, add it to _activeEdgesIneq
      _activeEdgesIneq.push_back(static_cast<OptimizableGraph::Edge *>(
          *it)); // Dereference the iterator to get the element
    if (_edgeEqSet.find(static_cast<OptimizableGraph::Edge *>(*it)) !=
        _edgeEqSet.end())
      _activeEdgesEq.push_back(static_cast<OptimizableGraph::Edge *>(*it));

    //-----------------------------------------
  }

  _activeVertices.reserve(auxVertexSet.size());
  for (set<Vertex *>::iterator it = auxVertexSet.begin();
       it != auxVertexSet.end(); ++it)
    _activeVertices.push_back(*it);

  sortVectorContainers();
  bool indexMappingStatus = buildIndexMapping(_activeVertices);
  postIteration(-1);
  return indexMappingStatus;
}

// Implementation of `addEdgeIneq`
template <typename Derived>
template <int D, typename E, typename... VertexTypes>
bool SparseOptimizerForConstraints<Derived>::addEdgeIneq(
    BaseFixedSizedEdgeIneq<D, E, VertexTypes...> *e) {
  return static_cast<Derived *>(this)->addEdgeIneqImpl(e);
}

// Implementation of `addEdgeEq`
template <typename Derived>
template <int D, typename E, typename... VertexTypes>
bool SparseOptimizerForConstraints<Derived>::addEdgeEq(
    BaseFixedSizedEdgeEq<D, E, VertexTypes...> *e) {
  return static_cast<Derived *>(this)->addEdgeEqImpl(e);
}

// Equality edges get a quadratic-form callback that applies the optimizer's
// default dual regularization unless the edge has its own value.
template <typename Derived>
template <int D, typename E, typename... VertexTypes>
void SparseOptimizerForConstraints<Derived>::installEqualityQuadraticForm(
    BaseFixedSizedEdgeEq<D, E, VertexTypes...> *e) {
  e->setConstructQuadraticFormImpl(
      [this](BaseFixedSizedEdgeEq<D, E, VertexTypes...> &edge) {
        const double own = edge.dualRegularization();
        edge.constructQuadraticFormRegularized(own >= 0.0 ? own
                                                          : _eq_regularization);
      });
}

// Termination functions
template <typename Derived>
TerminationCriteria &
SparseOptimizerForConstraints<Derived>::terminationCriteria() {
  return _terminationCriteria;
}

template <typename Derived>
void SparseOptimizerForConstraints<Derived>::setConvergenceCriterion(
    int criterion) {
  /*
  the termination criterion is based on the Newton decrement, the gradient norm,
  or the update norm; 0: UpdateNorm, 1: NewtonDecrement, 2: GradientNorm,
  3: RelativeNewtonDecrement (|x^T b| < eps (1 + cost))
  */
  _terminationCriteria.setConvergenceCriterion(criterion);
};

template <typename Derived>
bool SparseOptimizerForConstraints<Derived>::verifyConvergence(double epsilon) {
  g2o::OptimizationAlgorithmWithHessian *algorithm =
      static_cast<g2o::OptimizationAlgorithmWithHessian *>(_algorithm);

  // Get the gradient vector b and the vector size
  const double *b = algorithm->solver().b();
  const double *update = algorithm->solver().x();
  size_t xSize = algorithm->solver().vectorSize();

  // Map the b vector to an Eigen vector
  Eigen::Map<const Eigen::VectorXd> bVec(b, xSize);
  Eigen::Map<const Eigen::VectorXd> updateVec(update, xSize);

  // the objective value is only evaluated for the relative criterion
  double cost = 0.0;
  if (terminationCriteria().needsCost()) {
    computeActiveErrors();
    cost = this->activeRobustChi2();
  }

  if (terminationCriteria().verifyConvergence(bVec, updateVec, epsilon, cost)) {
    return true;
  }
  return false;
}

template <typename Derived>
bool SparseOptimizerForConstraints<Derived>::verifyIneqFeasibility(
    EdgeContainer &edgeVec, double epsilon) {
  return terminationCriteria().verifyIneqFeasibility(edgeVec, epsilon);
}

template <typename Derived>
bool SparseOptimizerForConstraints<Derived>::verifyEqFeasibility(
    EdgeContainer &edgeVec, double epsilon) {
  return terminationCriteria().verifyEqFeasibility(edgeVec, epsilon);
}

template <typename Derived>
void SparseOptimizerForConstraints<Derived>::setWarmStartLagrangeMultiplierFlag(
    bool flag) {
  static_cast<Derived *>(this)->setMultiplierReuse(flag);
}

template <typename Derived>
void SparseOptimizerForConstraints<Derived>::setInnerIterationsMax(
    int maxInnerIterations) {
  if (maxInnerIterations <= 0)
    throw std::invalid_argument("inner iteration limit must be positive");
  _num_inner_iterations_max = maxInnerIterations;
}

template <typename Derived>
void SparseOptimizerForConstraints<Derived>::setLagrangeMultiplierInitial(
    double lagrangeMultiplierInitial) {
  _lagrange_multiplier_initial_eq = lagrangeMultiplierInitial;
}

template <typename Derived>
void SparseOptimizerForConstraints<Derived>::setEpsilonConstraint(
    double epsilon_feas) {
  if (!(epsilon_feas >= 0) || !std::isfinite(epsilon_feas))
    throw std::invalid_argument(
        "constraint tolerance must be finite and nonnegative");
  _terminationCriteria.epsilon_constraint = epsilon_feas;
}

template <typename Derived>
void SparseOptimizerForConstraints<Derived>::setEpsilonConvergence(
    double epsilon) {
  if (!(epsilon >= 0) || !std::isfinite(epsilon))
    throw std::invalid_argument(
        "convergence tolerance must be finite and nonnegative");
  _terminationCriteria.epsilon_convergence = epsilon;
}

template <typename Derived>
void SparseOptimizerForConstraints<Derived>::setAlphaBacktrackingVector(
    std::vector<double> alphaBacktracking) {
  if (alphaBacktracking.empty())
    throw std::invalid_argument("backtracking sequence must not be empty");
  double previous = 2.0;
  for (double value : alphaBacktracking) {
    if (!(value > 0.0 && value <= 1.0 && value < previous) ||
        !std::isfinite(value))
      throw std::invalid_argument("backtracking steps must be finite, in "
                                  "(0,1], and strictly decreasing");
    previous = value;
  }
  _alphaBacktracking = std::move(alphaBacktracking);
}

template <typename Derived>
void SparseOptimizerForConstraints<Derived>::setAlphaBacktracking(
    double alpha) {
  defineAlphaBacktracking(alpha);
  _alpha = alpha;
}

template <typename Derived>
void SparseOptimizerForConstraints<Derived>::setNumInnerIterationsMax(
    int numInnerIterationsMax) {
  setInnerIterationsMax(numInnerIterationsMax);
}

template <typename Derived>
bool SparseOptimizerForConstraints<
    Derived>::getWarmStartLagrangeMultiplierFlag() {
  return static_cast<Derived *>(this)->multiplierReuse();
}

template <typename Derived>
int SparseOptimizerForConstraints<Derived>::getInnerIterationsMax() {
  return _num_inner_iterations_max;
}

template <typename Derived>
double SparseOptimizerForConstraints<Derived>::getLagrangeMultiplierInitial() {
  return _lagrange_multiplier_initial_eq;
}

template <typename Derived>
double SparseOptimizerForConstraints<Derived>::alpha() {
  return _alpha;
}

template <typename Derived>
double SparseOptimizerForConstraints<Derived>::NumInnerIterationsMax() {
  return _num_inner_iterations_max;
}

template <typename Derived>
void SparseOptimizerForConstraints<Derived>::defineAlphaBacktracking(
    double alpha) {
  if (alpha == -1.0) {
    _alphaBacktracking = {
        1,     0.9,   0.8,   0.7,   0.6,   0.5,   0.4,   0.3,   0.25,  0.2,
        0.15,  0.1,   0.09,  0.08,  0.06,  0.05,  0.04,  0.03,  0.02,  0.01,
        0.008, 0.007, 0.006, 0.005, 0.004, 0.003, 0.002, 1e-3,  1e-4,  1e-5,
        1e-7,  1e-9,  1e-10, 1e-12, 1e-14, 1e-16, 1e-18, 1e-20, 1e-25, 1e-30};
    return;
  } // the default case

  if (alpha > 0 && alpha < 1) {
    const double steps = std::ceil(std::log(1e-30) / std::log(alpha));
    if (!std::isfinite(steps) || steps > 100000)
      throw std::invalid_argument(
          "backtracking ratio would require more than 100000 steps");
    std::vector<double> sequence{1.0};
    sequence.reserve(static_cast<std::size_t>(steps) + 2);
    while (sequence.back() > 1e-30) {
      sequence.push_back(sequence.back() * alpha);
    }
    _alphaBacktracking.swap(sequence);
  } else {
    throw std::invalid_argument(
        "backtracking ratio must be in (0,1), or -1 for the standard sequence");
  }
}

} // namespace cg2o
