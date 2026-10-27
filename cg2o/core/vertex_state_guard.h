// Exact, temporary trial-state rollback for arbitrary g2o retractions.
#ifndef CG2O_VERTEX_STATE_GUARD_H
#define CG2O_VERTEX_STATE_GUARD_H

#include <g2o/core/optimizable_graph.h>
#include <algorithm>

namespace cg2o::internal {
// No auxiliary heap storage. Typical constraint arities are small. A vertex
// repeated in an edge is backed up once, and fixed measurements are skipped.
class VertexStateGuard {
public:
  explicit VertexStateGuard(const g2o::HyperGraph::VertexContainer& vertices)
      : _vertices(vertices) {
    for (std::size_t i = 0; i < _vertices.size(); ++i)
      if (eligible(i)) asVertex(i)->push();
  }
  VertexStateGuard(const VertexStateGuard&) = delete;
  VertexStateGuard& operator=(const VertexStateGuard&) = delete;
  ~VertexStateGuard() { restore(); }
  void restore() {
    if (!_active) return;
    for (std::size_t i = _vertices.size(); i > 0; --i)
      if (eligible(i-1)) asVertex(i-1)->pop();
    _active = false;
  }
private:
  g2o::OptimizableGraph::Vertex* asVertex(std::size_t i) const {
    return static_cast<g2o::OptimizableGraph::Vertex*>(_vertices[i]);
  }
  bool eligible(std::size_t i) const {
    return _vertices[i] && !asVertex(i)->fixed() &&
           std::find(_vertices.begin(), _vertices.begin()+i, _vertices[i]) == _vertices.begin()+i;
  }
  const g2o::HyperGraph::VertexContainer& _vertices;
  bool _active = true;
};
} // namespace cg2o::internal
#endif
