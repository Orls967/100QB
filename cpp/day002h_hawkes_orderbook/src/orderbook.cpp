#include "hawkes_ob/orderbook.hpp"

#include <algorithm>
#include <bit>

namespace hob {

OrderBook::OrderBook(std::int64_t min_price, std::int64_t max_price, std::size_t max_orders)
    : min_price_(min_price),
      max_price_(max_price),
      best_bid_(min_price - 1),
      best_ask_(max_price + 1),
      bid_levels_(static_cast<std::size_t>(max_price - min_price + 1)),
      ask_levels_(static_cast<std::size_t>(max_price - min_price + 1)),
      bid_bm1_((bid_levels_.size() + 63) / 64, 0ULL),
      bid_bm2_((bid_bm1_.size() + 63) / 64, 0ULL),
      ask_bm1_((ask_levels_.size() + 63) / 64, 0ULL),
      ask_bm2_((ask_bm1_.size() + 63) / 64, 0ULL),
      pool_(max_orders),
      free_list_(max_orders),
      id_to_slot_(max_orders, OrderNode::kNull),
      id_to_active_pos_(max_orders, OrderNode::kNull) {
  active_orders_.reserve(max_orders);
  for (std::size_t i = 0; i < max_orders; ++i) {
    free_list_[i] = static_cast<std::uint32_t>(max_orders - 1 - i);
  }
}

void OrderBook::reset() noexcept {
  best_bid_ = min_price_ - 1;
  best_ask_ = max_price_ + 1;
  for (auto& lvl : bid_levels_) {
    lvl = PriceLevel{};
  }
  for (auto& lvl : ask_levels_) {
    lvl = PriceLevel{};
  }
  std::fill(bid_bm1_.begin(), bid_bm1_.end(), 0ULL);
  std::fill(bid_bm2_.begin(), bid_bm2_.end(), 0ULL);
  std::fill(ask_bm1_.begin(), ask_bm1_.end(), 0ULL);
  std::fill(ask_bm2_.begin(), ask_bm2_.end(), 0ULL);
  for (auto& n : pool_) {
    n = OrderNode{};
  }
  const std::size_t n = pool_.size();
  free_list_.resize(n);
  for (std::size_t i = 0; i < n; ++i) {
    free_list_[i] = static_cast<std::uint32_t>(n - 1 - i);
  }
  std::fill(id_to_slot_.begin(), id_to_slot_.end(), OrderNode::kNull);
  std::fill(id_to_active_pos_.begin(), id_to_active_pos_.end(), OrderNode::kNull);
  active_orders_.clear();
}

void OrderBook::set_level_bit(Side side, std::size_t idx) noexcept {
  const std::size_t w1 = idx >> 6;
  const std::size_t b1 = idx & 63U;
  const std::size_t w2 = w1 >> 6;
  const std::size_t b2 = w1 & 63U;
  if (side == Side::Buy) {
    bid_bm1_[w1] |= (1ULL << b1);
    bid_bm2_[w2] |= (1ULL << b2);
  } else {
    ask_bm1_[w1] |= (1ULL << b1);
    ask_bm2_[w2] |= (1ULL << b2);
  }
}

void OrderBook::clear_level_bit(Side side, std::size_t idx) noexcept {
  const std::size_t w1 = idx >> 6;
  const std::size_t b1 = idx & 63U;
  const std::size_t w2 = w1 >> 6;
  const std::size_t b2 = w1 & 63U;
  if (side == Side::Buy) {
    bid_bm1_[w1] &= ~(1ULL << b1);
    if (bid_bm1_[w1] == 0ULL) {
      bid_bm2_[w2] &= ~(1ULL << b2);
    }
  } else {
    ask_bm1_[w1] &= ~(1ULL << b1);
    if (ask_bm1_[w1] == 0ULL) {
      ask_bm2_[w2] &= ~(1ULL << b2);
    }
  }
}

std::uint32_t OrderBook::alloc_slot() noexcept {
  const std::uint32_t slot = free_list_.back();
  free_list_.pop_back();
  return slot;
}

void OrderBook::free_slot(std::uint32_t slot) noexcept {
  const std::uint64_t id = pool_[slot].id;
  pool_[slot].active = false;
  pool_[slot].prev = OrderNode::kNull;
  pool_[slot].next = OrderNode::kNull;
  id_to_slot_[id] = OrderNode::kNull;

  const std::uint32_t pos = id_to_active_pos_[id];
  const std::uint64_t last_id = active_orders_.back();
  active_orders_[pos] = last_id;
  id_to_active_pos_[last_id] = pos;
  active_orders_.pop_back();
  id_to_active_pos_[id] = OrderNode::kNull;

  free_list_.push_back(slot);
}

void OrderBook::link_tail(PriceLevel& lvl, std::uint32_t slot) noexcept {
  pool_[slot].prev = lvl.tail;
  pool_[slot].next = OrderNode::kNull;
  if (lvl.tail != OrderNode::kNull) {
    pool_[lvl.tail].next = slot;
  } else {
    lvl.head = slot;
  }
  lvl.tail = slot;
  lvl.total_qty += pool_[slot].qty;
  ++lvl.count;
}

void OrderBook::unlink_node(PriceLevel& lvl, std::uint32_t slot) noexcept {
  const std::uint32_t p = pool_[slot].prev;
  const std::uint32_t n = pool_[slot].next;
  if (p != OrderNode::kNull) {
    pool_[p].next = n;
  } else {
    lvl.head = n;
  }
  if (n != OrderNode::kNull) {
    pool_[n].prev = p;
  } else {
    lvl.tail = p;
  }
  lvl.total_qty -= pool_[slot].qty;
  --lvl.count;
}

void OrderBook::advance_best_bid() noexcept {
  for (std::size_t w2 = bid_bm2_.size(); w2-- > 0;) {
    const std::uint64_t v2 = bid_bm2_[w2];
    if (v2 == 0ULL) {
      continue;
    }
    const std::size_t b2 = 63U - static_cast<std::size_t>(std::countl_zero(v2));
    const std::size_t w1 = (w2 << 6) | b2;
    const std::uint64_t v1 = bid_bm1_[w1];
    const std::size_t b1 = 63U - static_cast<std::size_t>(std::countl_zero(v1));
    best_bid_ = min_price_ + static_cast<std::int64_t>((w1 << 6) | b1);
    return;
  }
  best_bid_ = min_price_ - 1;
}

void OrderBook::advance_best_ask() noexcept {
  for (std::size_t w2 = 0; w2 < ask_bm2_.size(); ++w2) {
    const std::uint64_t v2 = ask_bm2_[w2];
    if (v2 == 0ULL) {
      continue;
    }
    const std::size_t b2 = static_cast<std::size_t>(std::countr_zero(v2));
    const std::size_t w1 = (w2 << 6) | b2;
    const std::uint64_t v1 = ask_bm1_[w1];
    const std::size_t b1 = static_cast<std::size_t>(std::countr_zero(v1));
    best_ask_ = min_price_ + static_cast<std::int64_t>((w1 << 6) | b1);
    return;
  }
  best_ask_ = max_price_ + 1;
}

void OrderBook::submit_limit(std::uint64_t id, Side side, std::int64_t price, std::int64_t qty,
                             std::vector<Trade>& out_trades) noexcept {
  out_trades.clear();
  if (qty <= 0 || price < min_price_ || price > max_price_) {
    return;
  }

  std::int64_t rem = qty;
  if (side == Side::Buy) {
    while (rem > 0 && best_ask_ <= max_price_ && best_ask_ <= price) {
      const std::size_t p_i = price_idx(best_ask_);
      PriceLevel& lvl = ask_levels_[p_i];
      while (rem > 0 && lvl.head != OrderNode::kNull) {
        const std::uint32_t maker_slot = lvl.head;
        OrderNode& maker = pool_[maker_slot];
        const std::int64_t fill = std::min(rem, maker.qty);
        out_trades.push_back({maker.id, id, best_ask_, fill});
        rem -= fill;
        maker.qty -= fill;
        lvl.total_qty -= fill;
        if (maker.qty == 0) {
          unlink_node(lvl, maker_slot);
          free_slot(maker_slot);
        }
      }
      if (lvl.count == 0) {
        clear_level_bit(Side::Sell, p_i);
        advance_best_ask();
      }
    }
    if (rem > 0) {
      const std::uint32_t slot = alloc_slot();
      pool_[slot] = {id, price, rem, Side::Buy, OrderNode::kNull, OrderNode::kNull, true};
      id_to_slot_[id] = slot;
      id_to_active_pos_[id] = static_cast<std::uint32_t>(active_orders_.size());
      active_orders_.push_back(id);
      const std::size_t p_i = price_idx(price);
      if (bid_levels_[p_i].count == 0) {
        set_level_bit(Side::Buy, p_i);
      }
      link_tail(bid_levels_[p_i], slot);
      if (price > best_bid_) {
        best_bid_ = price;
      }
    }
  } else {
    while (rem > 0 && best_bid_ >= min_price_ && best_bid_ >= price) {
      const std::size_t p_i = price_idx(best_bid_);
      PriceLevel& lvl = bid_levels_[p_i];
      while (rem > 0 && lvl.head != OrderNode::kNull) {
        const std::uint32_t maker_slot = lvl.head;
        OrderNode& maker = pool_[maker_slot];
        const std::int64_t fill = std::min(rem, maker.qty);
        out_trades.push_back({maker.id, id, best_bid_, fill});
        rem -= fill;
        maker.qty -= fill;
        lvl.total_qty -= fill;
        if (maker.qty == 0) {
          unlink_node(lvl, maker_slot);
          free_slot(maker_slot);
        }
      }
      if (lvl.count == 0) {
        clear_level_bit(Side::Buy, p_i);
        advance_best_bid();
      }
    }
    if (rem > 0) {
      const std::uint32_t slot = alloc_slot();
      pool_[slot] = {id, price, rem, Side::Sell, OrderNode::kNull, OrderNode::kNull, true};
      id_to_slot_[id] = slot;
      id_to_active_pos_[id] = static_cast<std::uint32_t>(active_orders_.size());
      active_orders_.push_back(id);
      const std::size_t p_i = price_idx(price);
      if (ask_levels_[p_i].count == 0) {
        set_level_bit(Side::Sell, p_i);
      }
      link_tail(ask_levels_[p_i], slot);
      if (price < best_ask_) {
        best_ask_ = price;
      }
    }
  }
}

void OrderBook::submit_market(std::uint64_t id, Side side, std::int64_t qty,
                              std::vector<Trade>& out_trades) noexcept {
  out_trades.clear();
  if (qty <= 0) {
    return;
  }
  std::int64_t rem = qty;
  if (side == Side::Buy) {
    while (rem > 0 && best_ask_ <= max_price_) {
      const std::size_t p_i = price_idx(best_ask_);
      PriceLevel& lvl = ask_levels_[p_i];
      while (rem > 0 && lvl.head != OrderNode::kNull) {
        const std::uint32_t maker_slot = lvl.head;
        OrderNode& maker = pool_[maker_slot];
        const std::int64_t fill = std::min(rem, maker.qty);
        out_trades.push_back({maker.id, id, best_ask_, fill});
        rem -= fill;
        maker.qty -= fill;
        lvl.total_qty -= fill;
        if (maker.qty == 0) {
          unlink_node(lvl, maker_slot);
          free_slot(maker_slot);
        }
      }
      if (lvl.count == 0) {
        clear_level_bit(Side::Sell, p_i);
        advance_best_ask();
      }
    }
  } else {
    while (rem > 0 && best_bid_ >= min_price_) {
      const std::size_t p_i = price_idx(best_bid_);
      PriceLevel& lvl = bid_levels_[p_i];
      while (rem > 0 && lvl.head != OrderNode::kNull) {
        const std::uint32_t maker_slot = lvl.head;
        OrderNode& maker = pool_[maker_slot];
        const std::int64_t fill = std::min(rem, maker.qty);
        out_trades.push_back({maker.id, id, best_bid_, fill});
        rem -= fill;
        maker.qty -= fill;
        lvl.total_qty -= fill;
        if (maker.qty == 0) {
          unlink_node(lvl, maker_slot);
          free_slot(maker_slot);
        }
      }
      if (lvl.count == 0) {
        clear_level_bit(Side::Buy, p_i);
        advance_best_bid();
      }
    }
  }
}

bool OrderBook::cancel(std::uint64_t id) noexcept {
  if (id >= id_to_slot_.size()) {
    return false;
  }
  const std::uint32_t slot = id_to_slot_[id];
  if (slot == OrderNode::kNull) {
    return false;
  }
  const Side side = pool_[slot].side;
  const std::int64_t price = pool_[slot].price;
  const std::size_t p_i = price_idx(price);
  PriceLevel& lvl = (side == Side::Buy) ? bid_levels_[p_i] : ask_levels_[p_i];
  unlink_node(lvl, slot);
  free_slot(slot);
  if (lvl.count == 0) {
    clear_level_bit(side, p_i);
    if (side == Side::Buy && price == best_bid_) {
      advance_best_bid();
    } else if (side == Side::Sell && price == best_ask_) {
      advance_best_ask();
    }
  }
  return true;
}

bool OrderBook::modify(std::uint64_t id, std::int64_t new_qty) noexcept {
  if (new_qty <= 0) {
    return cancel(id);
  }
  if (id >= id_to_slot_.size()) {
    return false;
  }
  const std::uint32_t slot = id_to_slot_[id];
  if (slot == OrderNode::kNull) {
    return false;
  }
  OrderNode& node = pool_[slot];
  const std::size_t p_i = price_idx(node.price);
  PriceLevel& lvl = (node.side == Side::Buy) ? bid_levels_[p_i] : ask_levels_[p_i];
  if (new_qty <= node.qty) {
    lvl.total_qty -= (node.qty - new_qty);
    node.qty = new_qty;
  } else {
    unlink_node(lvl, slot);
    node.qty = new_qty;
    link_tail(lvl, slot);
  }
  return true;
}

std::vector<LevelSnapshot> OrderBook::snapshot() const {
  std::vector<LevelSnapshot> snap;
  for (std::int64_t p = best_bid_; p >= min_price_; --p) {
    const auto& lvl = bid_levels_[price_idx(p)];
    if (lvl.count > 0) {
      snap.push_back({Side::Buy, p, lvl.total_qty, lvl.count});
    }
  }
  for (std::int64_t p = best_ask_; p <= max_price_; ++p) {
    const auto& lvl = ask_levels_[price_idx(p)];
    if (lvl.count > 0) {
      snap.push_back({Side::Sell, p, lvl.total_qty, lvl.count});
    }
  }
  return snap;
}

bool OrderBook::check_invariants() const noexcept {
  if (best_bid_ >= min_price_ && best_ask_ <= max_price_ && best_bid_ >= best_ask_) {
    return false;
  }
  std::size_t total_orders = 0;
  for (std::int64_t p = min_price_; p <= max_price_; ++p) {
    const auto& bl = bid_levels_[price_idx(p)];
    if (bl.count > 0 && p > best_bid_) {
      return false;
    }
    std::int64_t sum_q = 0;
    std::size_t cnt = 0;
    for (std::uint32_t s = bl.head; s != OrderNode::kNull; s = pool_[s].next) {
      if (!pool_[s].active || pool_[s].qty <= 0 || id_to_slot_[pool_[s].id] != s) {
        return false;
      }
      sum_q += pool_[s].qty;
      ++cnt;
    }
    if (sum_q != bl.total_qty || cnt != bl.count) {
      return false;
    }
    total_orders += cnt;

    const auto& al = ask_levels_[price_idx(p)];
    if (al.count > 0 && p < best_ask_) {
      return false;
    }
    sum_q = 0;
    cnt = 0;
    for (std::uint32_t s = al.head; s != OrderNode::kNull; s = pool_[s].next) {
      if (!pool_[s].active || pool_[s].qty <= 0 || id_to_slot_[pool_[s].id] != s) {
        return false;
      }
      sum_q += pool_[s].qty;
      ++cnt;
    }
    if (sum_q != al.total_qty || cnt != al.count) {
      return false;
    }
    total_orders += cnt;
  }
  return total_orders == active_orders_.size();
}

} // namespace hob
