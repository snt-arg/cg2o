// Building blocks shared by the start-point policies of the constrained
// optimizers: solve outcome, the copies of vertex estimates, stage shifts.
#ifndef CG2O_START_POINT_COMMON_H
#define CG2O_START_POINT_COMMON_H

#include <algorithm>
#include <cstddef>
#include <limits>
#include <memory>
#include <optimizable_graph.h>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace cg2o {

// The rows of one inequality edge, independent of its dimension D: slack s
// and multiplier lambda live in the edge; g(x) is edge->errorData() after
// edge->computeError().
struct IneqBlock {
  g2o::OptimizableGraph::Edge *edge = nullptr;
  int dimension = 0;
  double *slack = nullptr;
  double *multiplier = nullptr;
};

//! what a start-point policy reports to optimize() before the iterations
struct StartResult {
  bool warm = false;         //!< started from cached or supplied values
  bool keepRollback = false; //!< the starting point was stored for rollback
  //! max(|h|, g+) at the start, measured for the divergence test
  double startViolation = std::numeric_limits<double>::quiet_NaN();
};

//! the result of one run of the iterations
struct SolveOutcome {
  int iterations = 0;        //!< -1: the algorithm failed to initialize
  bool solverFailed = false; //!< the linear solver failed in an iteration
};

//! what optimize() does after one run of the iterations
enum class FinishAction {
  Accept,    //!< the iterate is the result
  RetryCold, //!< the starting point was put back: solve again, cold
  Reject,    //!< the solve broke down; the cache is cleared
  Error      //!< the algorithm failed to initialize (optimize() returns -1)
};

//! stage k takes stage min(k + n, last): the last n stages keep the last one
inline std::size_t shiftSource(std::size_t k, std::size_t n, std::size_t size) {
  return std::min(k + n, size - 1);
}

// Copies of vertex estimates for the start points of a solver: the cached
// point of the next solve and the rollback point of the current one. One
// buffer per vertex and copy is allocated when the vertex is registered and
// reused afterwards. Estimates are copied as estimates (manifold-safe),
// never as increments, and never through g2o's push/pop backup stacks.
class VertexEstimateCache {
public:
  enum Copy { Cached = 0, Rollback = 1 };

  //! typed copies through estimate()/setEstimate(); replaces earlier ones
  template <typename T> void registerVertex(T *v) {
    replace(v, std::make_unique<TypedSlot<T>>(v));
  }
  void forgetVertex(const g2o::OptimizableGraph::Vertex *v) {
    auto it = _slots.find(v);
    if (it == _slots.end())
      return;
    retire(std::move(it->second));
    _slots.erase(it);
  }
  void clear() {
    discard(Cached);
    discard(Rollback);
    _slots.clear();
    _retired.clear();
  }

  //! copies the non-fixed active vertices: the equality multipliers always,
  //! the other vertices only with `primal`. A vertex without typed
  //! registration needs the estimate-data interface.
  void store(Copy which,
             const std::vector<g2o::OptimizableGraph::Vertex *> &active,
             const std::vector<g2o::OptimizableGraph::Vertex *> &multipliers,
             bool primal) {
    discard(which); // at most one entry per vertex and copy
    _multipliers.assign(multipliers.begin(), multipliers.end());
    std::sort(_multipliers.begin(), _multipliers.end());
    auto &entries = _entries[which];
    entries.reserve(active.size());
    for (auto *v : active) {
      if (v->fixed())
        continue;
      const bool multiplier =
          std::binary_search(_multipliers.begin(), _multipliers.end(), v);
      if (!multiplier && !primal)
        continue; // primal not restored, so not stored
      auto it = _slots.find(v);
      if (it == _slots.end())
        it = _slots.emplace(v, std::make_unique<DataSlot>(v)).first;
      if (!it->second->save(which))
        throw std::invalid_argument("warm start needs typed addVertex() or "
                                    "registerWarmStartVertex<T>() for vertex " +
                                    std::to_string(v->id()));
      entries.push_back({v, v->id(), multiplier, it->second.get()});
    }
  }

  //! puts the copies back into the vertices that are still in the graph and
  //! active; fixed vertices keep their current (measured) values
  void restore(Copy which, g2o::OptimizableGraph &graph, bool primal) {
    for (auto &entry : _entries[which]) {
      auto *v = graph.vertex(entry.id);
      if (v != entry.vertex || v->hessianIndex() < 0)
        continue;
      if (!v->fixed() && (entry.multiplier || primal)) {
        if (!entry.slot->restore(which))
          throw std::runtime_error(
              "failed to restore warm-start vertex estimate");
      }
    }
    discard(which);
  }

  void discard(Copy which) {
    _entries[which].clear();
    if (_entries[Cached].empty() && _entries[Rollback].empty())
      _retired.clear();
  }

private:
  struct Slot {
    virtual ~Slot() = default;
    virtual bool save(int which) = 0;
    virtual bool restore(int which) = 0;
  };
  template <typename T> struct TypedSlot : Slot {
    using Estimate = std::decay_t<decltype(std::declval<T &>().estimate())>;
    explicit TypedSlot(T *vertex) : v(vertex) {}
    bool save(int which) override {
      saved[which] = v->estimate();
      return true;
    }
    bool restore(int which) override {
      v->setEstimate(saved[which]);
      return true;
    }
    T *v;
    Estimate saved[2];
  };
  struct DataSlot : Slot {
    explicit DataSlot(g2o::OptimizableGraph::Vertex *vertex) : v(vertex) {}
    bool save(int which) override {
      return v->estimateDimension() > 0 && v->getEstimateData(saved[which]);
    }
    bool restore(int which) override {
      return v->setEstimateData(saved[which].data());
    }
    g2o::OptimizableGraph::Vertex *v;
    std::vector<double> saved[2];
  };
  struct Entry {
    g2o::OptimizableGraph::Vertex *vertex;
    int id;
    bool multiplier; // equality multiplier vertex (nu), not primal
    Slot *slot;
  };
  void replace(const g2o::OptimizableGraph::Vertex *v,
               std::unique_ptr<Slot> slot) {
    auto &current = _slots[v];
    if (current)
      retire(std::move(current));
    current = std::move(slot);
  }
  // a stored copy may still refer to a replaced slot: keep it until both
  // copies are discarded
  void retire(std::unique_ptr<Slot> slot) {
    if (!_entries[Cached].empty() || !_entries[Rollback].empty())
      _retired.push_back(std::move(slot));
  }

  std::unordered_map<const g2o::OptimizableGraph::Vertex *,
                     std::unique_ptr<Slot>>
      _slots;
  std::vector<std::unique_ptr<Slot>> _retired;
  std::vector<Entry> _entries[2];
  std::vector<g2o::OptimizableGraph::Vertex *> _multipliers; // sorted
};

// Disjoint stage sequences registered for the shift of a receding horizon.
class ShiftMembers {
public:
  //! rejects null, repeated and already registered entries
  template <typename T> void check(const std::vector<T *> &stages) const {
    for (auto *p : stages) {
      if (p == nullptr)
        throw std::invalid_argument(
            "[Error] addWarmStartShiftSequence: null entry in the sequence");
      if (_members.count(p))
        throw std::invalid_argument("warm-start shift sequences overlap");
    }
    std::unordered_set<T *> unique(stages.begin(), stages.end());
    if (unique.size() != stages.size())
      throw std::invalid_argument(
          "duplicate entry in warm-start shift sequence");
  }
  template <typename T> void add(const std::vector<T *> &stages) {
    _members.insert(stages.begin(), stages.end());
  }
  void clear() { _members.clear(); }

private:
  std::unordered_set<const void *> _members;
};

} // namespace cg2o

#endif // CG2O_START_POINT_COMMON_H
