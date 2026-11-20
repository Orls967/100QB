#pragma once

#include <cstdint>
#include <limits>
#include <vector>

namespace hob {

enum class Side : std::uint8_t { Buy = 0, Sell = 1 };

struct Trade {
  std::uint64_t maker_id{0};
  std::uint64_t taker_id{0};
  std::int64_t price{0};
  std::int64_t qty{0};

  friend bool operator==(const Trade&, const Trade&) = default;
};

struct LevelSnapshot {
  Side side{Side::Buy};
  std::int64_t price{0};
  std::int64_t qty{0};
  std::size_t count{0};

  friend bool operator==(const LevelSnapshot&, const LevelSnapshot&) = default;
};

struct OrderNode {
  static constexpr std::uint32_t kNull = std::numeric_limits<std::uint32_t>::max();
  std::uint64_t id{0};
  std::int64_t price{0};
  std::int64_t qty{0};
  Side side{Side::Buy};
  std::uint32_t prev{kNull};
  std::uint32_t next{kNull};
  bool active{false};
};

struct PriceLevel {
  std::uint32_t head{OrderNode::kNull};
  std::uint32_t tail{OrderNode::kNull};
  std::int64_t total_qty{0};
  std::size_t count{0};
};

/// High-performance price-time priority Limit Order Book with pool allocator,
/// intrusive FIFO queues, O(1) id lookup, and zero heap allocations in hot path.
class OrderBook {
public:
  OrderBook(std::int64_t min_price, std::int64_t max_price, std::size_t max_orders);

  void reset() noexcept;

  void submit_limit(std::uint64_t id, Side side, std::int64_t price, std::int64_t qty,
                    std::vector<Trade>& out_trades) noexcept;

  void submit_market(std::uint64_t id, Side side, std::int64_t qty,
                     std::vector<Trade>& out_trades) noexcept;

  bool cancel(std::uint64_t id) noexcept;

  bool modify(std::uint64_t id, std::int64_t new_qty) noexcept;

  [[nodiscard]] std::int64_t best_bid() const noexcept { return best_bid_; }
  [[nodiscard]] std::int64_t best_ask() const noexcept { return best_ask_; }
  [[nodiscard]] std::size_t active_order_count() const noexcept { return active_orders_.size(); }
  [[nodiscard]] const std::vector<std::uint64_t>& active_ids() const noexcept {
    return active_orders_;
  }

  [[nodiscard]] std::vector<LevelSnapshot> snapshot() const;
  [[nodiscard]] bool check_invariants() const noexcept;

private:
  [[nodiscard]] std::size_t price_idx(std::int64_t price) const noexcept {
    return static_cast<std::size_t>(price - min_price_);
  }

  std::uint32_t alloc_slot() noexcept;
  void free_slot(std::uint32_t slot) noexcept;
  void link_tail(PriceLevel& lvl, std::uint32_t slot) noexcept;
  void unlink_node(PriceLevel& lvl, std::uint32_t slot) noexcept;
  void advance_best_bid() noexcept;
  void advance_best_ask() noexcept;

  void set_level_bit(Side side, std::size_t idx) noexcept;
  void clear_level_bit(Side side, std::size_t idx) noexcept;

  std::int64_t min_price_;
  std::int64_t max_price_;
  std::int64_t best_bid_;
  std::int64_t best_ask_;

  std::vector<PriceLevel> bid_levels_;
  std::vector<PriceLevel> ask_levels_;
  std::vector<std::uint64_t> bid_bm1_;
  std::vector<std::uint64_t> bid_bm2_;
  std::vector<std::uint64_t> ask_bm1_;
  std::vector<std::uint64_t> ask_bm2_;
  std::vector<OrderNode> pool_;
  std::vector<std::uint32_t> free_list_;
  std::vector<std::uint32_t> id_to_slot_;
  std::vector<std::uint32_t> id_to_active_pos_;
  std::vector<std::uint64_t> active_orders_;
};

} // namespace hob
