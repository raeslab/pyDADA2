// parallel_for replacing RcppParallel::parallelFor in the ported DADA2 core (spec §5.3).
//
// Work is split into contiguous chunks of `grain` indices that are claimed dynamically by
// `n_threads` std::threads. The body receives [begin, end) ranges. Results are deterministic
// as long as each index writes only to its own output slot (as in all upstream Workers).
// Exceptions thrown by the body are captured and the first one is rethrown in the caller.

#ifndef DADA2_PARALLEL_H
#define DADA2_PARALLEL_H

#include <algorithm>
#include <atomic>
#include <cstddef>
#include <exception>
#include <mutex>
#include <thread>
#include <vector>

namespace dada2 {

template <typename F>
void parallel_for(std::size_t begin, std::size_t end, F &&body, std::size_t grain, int n_threads) {
  if (end <= begin) return;
  if (grain < 1) grain = 1;
  const std::size_t n = end - begin;
  const std::size_t nchunks = (n + grain - 1) / grain;
  std::size_t nthreads = n_threads < 1 ? 1 : static_cast<std::size_t>(n_threads);
  nthreads = std::min(nthreads, nchunks);
  if (nthreads <= 1) {
    body(begin, end);
    return;
  }

  std::atomic<std::size_t> next{0};
  std::atomic<bool> failed{false};
  std::exception_ptr error;
  std::mutex error_mutex;

  auto worker = [&]() {
    for (;;) {
      if (failed.load(std::memory_order_relaxed)) return;
      const std::size_t c = next.fetch_add(1);
      if (c >= nchunks) return;
      const std::size_t b = begin + c * grain;
      const std::size_t e = std::min(end, b + grain);
      try {
        body(b, e);
      } catch (...) {
        std::lock_guard<std::mutex> lock(error_mutex);
        if (!error) error = std::current_exception();
        failed.store(true);
        return;
      }
    }
  };

  std::vector<std::thread> threads;
  threads.reserve(nthreads - 1);
  for (std::size_t t = 0; t + 1 < nthreads; ++t) threads.emplace_back(worker);
  worker();
  for (auto &t : threads) t.join();
  if (error) std::rethrow_exception(error);
}

}  // namespace dada2

#endif  // DADA2_PARALLEL_H
