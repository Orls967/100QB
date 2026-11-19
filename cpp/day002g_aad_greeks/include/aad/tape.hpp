#pragma once

#include <cstddef>
#include <limits>
#include <vector>

namespace aad {

struct TapeNode {
  double value{0.0};
  double adjoint{0.0};
  double d0{0.0};
  double d1{0.0};
  std::size_t p0{kNullIdx};
  std::size_t p1{kNullIdx};

  static constexpr std::size_t kNullIdx = std::numeric_limits<std::size_t>::max();
};

class Tape {
public:
  Tape() { nodes_.reserve(4096); }

  void reserve(std::size_t cap) { nodes_.reserve(cap); }

  void clear() noexcept {
    nodes_.clear();
    mark_ = 0;
  }

  void mark() noexcept { mark_ = nodes_.size(); }

  void rewind_to_mark() noexcept { nodes_.resize(mark_); }

  [[nodiscard]] std::size_t mark_index() const noexcept { return mark_; }

  [[nodiscard]] std::size_t size() const noexcept { return nodes_.size(); }

  [[nodiscard]] std::size_t capacity_bytes() const noexcept {
    return nodes_.capacity() * sizeof(TapeNode);
  }

  std::size_t record_leaf(double val) {
    const std::size_t idx = nodes_.size();
    nodes_.push_back(TapeNode{val, 0.0, 0.0, 0.0, TapeNode::kNullIdx, TapeNode::kNullIdx});
    return idx;
  }

  std::size_t record_unary(double val, std::size_t p0, double d0) {
    const std::size_t idx = nodes_.size();
    nodes_.push_back(TapeNode{val, 0.0, d0, 0.0, p0, TapeNode::kNullIdx});
    return idx;
  }

  std::size_t record_binary(double val, std::size_t p0, double d0, std::size_t p1, double d1) {
    const std::size_t idx = nodes_.size();
    nodes_.push_back(TapeNode{val, 0.0, d0, d1, p0, p1});
    return idx;
  }

  [[nodiscard]] double value(std::size_t idx) const noexcept { return nodes_[idx].value; }

  [[nodiscard]] double adjoint(std::size_t idx) const noexcept { return nodes_[idx].adjoint; }

  void zero_all_adjoints() noexcept {
    for (auto& node : nodes_) {
      node.adjoint = 0.0;
    }
  }

  void backprop(std::size_t root_idx) noexcept {
    zero_all_adjoints();
    nodes_[root_idx].adjoint = 1.0;
    for (std::size_t i = root_idx + 1; i-- > 0;) {
      const double adj = nodes_[i].adjoint;
      if (adj == 0.0) {
        continue;
      }
      const std::size_t p0 = nodes_[i].p0;
      if (p0 != TapeNode::kNullIdx) {
        nodes_[p0].adjoint += adj * nodes_[i].d0;
        const std::size_t p1 = nodes_[i].p1;
        if (p1 != TapeNode::kNullIdx) {
          nodes_[p1].adjoint += adj * nodes_[i].d1;
        }
      }
    }
  }

  /// Reverse sweep only from `root_idx` down to `mark_` (accumulates into `[0, mark_)`).
  void backprop_to_mark(std::size_t root_idx, double seed_adj = 1.0) noexcept {
    for (std::size_t i = mark_; i < nodes_.size(); ++i) {
      nodes_[i].adjoint = 0.0;
    }
    nodes_[root_idx].adjoint = seed_adj;
    for (std::size_t i = root_idx + 1; i-- > mark_;) {
      const double adj = nodes_[i].adjoint;
      if (adj == 0.0) {
        continue;
      }
      const std::size_t p0 = nodes_[i].p0;
      if (p0 != TapeNode::kNullIdx) {
        nodes_[p0].adjoint += adj * nodes_[i].d0;
        const std::size_t p1 = nodes_[i].p1;
        if (p1 != TapeNode::kNullIdx) {
          nodes_[p1].adjoint += adj * nodes_[i].d1;
        }
      }
    }
  }

  /// Final reverse sweep from `mark_ - 1` down to `0` after all paths accumulated into `[0, mark_)`.
  void backprop_from_mark() noexcept {
    if (mark_ == 0) {
      return;
    }
    for (std::size_t i = mark_; i-- > 0;) {
      const double adj = nodes_[i].adjoint;
      if (adj == 0.0) {
        continue;
      }
      const std::size_t p0 = nodes_[i].p0;
      if (p0 != TapeNode::kNullIdx) {
        nodes_[p0].adjoint += adj * nodes_[i].d0;
        const std::size_t p1 = nodes_[i].p1;
        if (p1 != TapeNode::kNullIdx) {
          nodes_[p1].adjoint += adj * nodes_[i].d1;
        }
      }
    }
  }

private:
  std::vector<TapeNode> nodes_;
  std::size_t mark_{0};
};

inline Tape*& active_tape_slot() noexcept {
  thread_local Tape default_tape;
  thread_local Tape* current = &default_tape;
  return current;
}

inline Tape& active_tape() noexcept { return *active_tape_slot(); }

class TapeGuard {
public:
  explicit TapeGuard(Tape& tape) noexcept : prev_(active_tape_slot()) { active_tape_slot() = &tape; }
  ~TapeGuard() { active_tape_slot() = prev_; }
  TapeGuard(const TapeGuard&) = delete;
  TapeGuard& operator=(const TapeGuard&) = delete;

private:
  Tape* prev_;
};

} // namespace aad
