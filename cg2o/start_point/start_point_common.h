// Building blocks shared by the start points of the constrained solvers: the
// selected start method, the result of a run of the iterations, copies of
// vertex estimates and the stage sequences of a receding-horizon shift.
#ifndef CG2O_START_POINT_COMMON_H
#define CG2O_START_POINT_COMMON_H

#include <algorithm>
#include <concepts>
#include <cstddef>
#include <limits>
#include <memory>
#include <optimizable_graph.h>
#include <optional>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <typeindex>
#include <typeinfo>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

namespace cg2o {

//! the result of one run of the iterations
struct SolveOutcome {
  int iterations = 0;        //!< -1: the algorithm failed to initialize
  bool solverFailed = false; //!< the linear solver failed in an iteration
};

// The start method of a solver, one of the strategies of its family (Base is
// the family's strategy interface). A strategy S has a Settings type and a
// settings member. A method is selected between solves; the iterations only
// call the active one.
template <typename Base> class StrategySlot {
public:
  explicit StrategySlot(std::unique_ptr<Base> initial)
      : _active(std::move(initial)) {}

  Base &active() { return *_active; }
  const Base &active() const { return *_active; }
  std::type_index type() const {
    const Base &active = *_active;
    return typeid(active);
  }
  template <typename S> bool holds() const { return type() == typeid(S); }

  //! makes S the active method; an active S keeps its settings
  template <typename S>
    requires std::derived_from<S, Base>
  typename S::Settings &use() {
    if (!holds<S>())
      _active = std::make_unique<S>();
    return static_cast<S &>(*_active).settings;
  }
  //! settings of the active method S; std::logic_error if S is not active
  template <typename S>
    requires std::derived_from<S, Base>
  typename S::Settings &settings() {
    if (!holds<S>())
      throw std::logic_error(
          "settings requested for a start method that is not selected");
    return static_cast<S &>(*_active).settings;
  }
  template <typename S>
    requires std::derived_from<S, Base>
  const typename S::Settings &settings() const {
    if (!holds<S>())
      throw std::logic_error(
          "settings requested for a start method that is not selected");
    return static_cast<const S &>(*_active).settings;
  }

private:
  std::unique_ptr<Base> _active;
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
        throw std::invalid_argument(
            "the start point needs typed addVertex() or "
            "startPoint().registerVertex<T>() for vertex " +
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
          throw std::runtime_error("failed to restore a vertex estimate");
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
      if (!saved[which])
        return false;
      v->setEstimate(*saved[which]);
      return true;
    }
    T *v;
    // optional: the estimate type need not be default-constructible
    std::optional<Estimate> saved[2];
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
            "addShiftSequence: null entry in the sequence");
      if (_members.count(p))
        throw std::invalid_argument("shift sequences overlap");
    }
    std::unordered_set<T *> unique(stages.begin(), stages.end());
    if (unique.size() != stages.size())
      throw std::invalid_argument("duplicate entry in a shift sequence");
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
