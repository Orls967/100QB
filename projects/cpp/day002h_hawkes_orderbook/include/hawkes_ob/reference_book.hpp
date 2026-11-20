#pragma once

#include <algorithm>
#include <cstdint>
#include <map>
#include <vector>

#include "hawkes_ob/orderbook.hpp"

namespace hob {

/// Naive O(n) reference Limit Order Book using sorted vectors for exact equivalence testing.
class ReferenceBook {
public:
  struct Entry {
    std::uint64_t id{0};
    Side side{Side::Buy};
    std::int64_t price{0};
    std::int64_t qty{0};
  };

  void submit_limit(std::uint64_t id, Side side, std::int64_t price, std::int64_t qty,
                    std::vector<Trade>& out_trades) {
    out_trades.clear();
    if (qty <= 0) {
      return;
    }
    std::int64_t rem = qty;
    if (side == Side::Buy) {
      while (rem > 0 && !asks_.empty() && asks_.front().price <= price) {
        auto& maker = asks_.front();
        const std::int64_t fill = std::min(rem, maker.qty);
        out_trades.push_back({maker.id, id, maker.price, fill});
        rem -= fill;
        maker.qty -= fill;
        if (maker.qty == 0) {
          asks_.erase(asks_.begin());
        }
      }
      if (rem > 0) {
        auto it = std::upper_bound(bids_.begin(), bids_.end(), price,
                                   [](std::int64_t p, const Entry& e) { return p > e.price; });
        bids_.insert(it, {id, Side::Buy, price, rem});
      }
    } else {
      while (rem > 0 && !bids_.empty() && bids_.front().price >= price) {
        auto& maker = bids_.front();
        const std::int64_t fill = std::min(rem, maker.qty);
        out_trades.push_back({maker.id, id, maker.price, fill});
        rem -= fill;
        maker.qty -= fill;
        if (maker.qty == 0) {
          bids_.erase(bids_.begin());
        }
      }
      if (rem > 0) {
        auto it = std::upper_bound(asks_.begin(), asks_.end(), price,
                                   [](std::int64_t p, const Entry& e) { return p < e.price; });
        asks_.insert(it, {id, Side::Sell, price, rem});
      }
    }
  }

  void submit_market(std::uint64_t id, Side side, std::int64_t qty,
                     std::vector<Trade>& out_trades) {
    out_trades.clear();
    if (qty <= 0) {
      return;
    }
    std::int64_t rem = qty;
    if (side == Side::Buy) {
      while (rem > 0 && !asks_.empty()) {
        auto& maker = asks_.front();
        const std::int64_t fill = std::min(rem, maker.qty);
        out_trades.push_back({maker.id, id, maker.price, fill});
        rem -= fill;
        maker.qty -= fill;
        if (maker.qty == 0) {
          asks_.erase(asks_.begin());
        }
      }
    } else {
      while (rem > 0 && !bids_.empty()) {
        auto& maker = bids_.front();
        const std::int64_t fill = std::min(rem, maker.qty);
        out_trades.push_back({maker.id, id, maker.price, fill});
        rem -= fill;
        maker.qty -= fill;
        if (maker.qty == 0) {
          bids_.erase(bids_.begin());
        }
      }
    }
  }

  bool cancel(std::uint64_t id) {
    for (auto it = bids_.begin(); it != bids_.end(); ++it) {
      if (it->id == id) {
        bids_.erase(it);
        return true;
      }
    }
    for (auto it = asks_.begin(); it != asks_.end(); ++it) {
      if (it->id == id) {
        asks_.erase(it);
        return true;
      }
    }
    return false;
  }

  bool modify(std::uint64_t id, std::int64_t new_qty) {
    if (new_qty <= 0) {
      return cancel(id);
    }
    auto mod_vec = [&](std::vector<Entry>& vec) -> bool {
      for (auto it = vec.begin(); it != vec.end(); ++it) {
        if (it->id == id) {
          if (new_qty <= it->qty) {
            it->qty = new_qty;
          } else {
            Entry e = *it;
            e.qty = new_qty;
            vec.erase(it);
            if (e.side == Side::Buy) {
              auto ins = std::upper_bound(vec.begin(), vec.end(), e.price,
                                          [](std::int64_t p, const Entry& x) { return p > x.price; });
              vec.insert(ins, e);
            } else {
              auto ins = std::upper_bound(vec.begin(), vec.end(), e.price,
                                          [](std::int64_t p, const Entry& x) { return p < x.price; });
              vec.insert(ins, e);
            }
          }
          return true;
        }
      }
      return false;
    };
    return mod_vec(bids_) || mod_vec(asks_);
  }

  [[nodiscard]] std::vector<LevelSnapshot> snapshot() const {
    std::vector<LevelSnapshot> snap;
    for (const auto& e : bids_) {
      if (!snap.empty() && snap.back().side == Side::Buy && snap.back().price == e.price) {
        snap.back().qty += e.qty;
        ++snap.back().count;
      } else {
        snap.push_back({Side::Buy, e.price, e.qty, 1});
      }
    }
    for (const auto& e : asks_) {
      if (!snap.empty() && snap.back().side == Side::Sell && snap.back().price == e.price) {
        snap.back().qty += e.qty;
        ++snap.back().count;
      } else {
        snap.push_back({Side::Sell, e.price, e.qty, 1});
      }
    }
    return snap;
  }

private:
  std::vector<Entry> bids_;
  std::vector<Entry> asks_;
};

} // namespace hob
