#pragma once

#include <cstddef>
#include <thread>
#include <vector>

namespace lbm {

// Deterministic static range partition [0, n) across `num_threads` std::threads.
// When num_threads <= 1, executes synchronously on caller thread.
template <typename Func>
inline void parallel_for_range(std::size_t n, int num_threads, Func&& fn) {
    if (n == 0) return;
    int t = num_threads < 1 ? 1 : num_threads;
    if (static_cast<std::size_t>(t) > n) t = static_cast<int>(n);
    if (t == 1) {
        fn(std::size_t{0}, n, 0);
        return;
    }

    std::vector<std::thread> workers;
    workers.reserve(static_cast<std::size_t>(t - 1));
    const std::size_t chunk = (n + static_cast<std::size_t>(t) - 1) / static_cast<std::size_t>(t);
    for (int tid = 1; tid < t; ++tid) {
        const std::size_t beg = static_cast<std::size_t>(tid) * chunk;
        const std::size_t end = (beg + chunk < n) ? (beg + chunk) : n;
        workers.emplace_back([&, beg, end, tid]() { fn(beg, end, tid); });
    }
    const std::size_t end0 = (chunk < n) ? chunk : n;
    fn(std::size_t{0}, end0, 0);
    for (auto& w : workers) w.join();
}

}  // namespace lbm
