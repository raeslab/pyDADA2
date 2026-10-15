// parallel_for replacing RcppParallel::parallelFor in the ported DADA2 core (spec §5.3).
//
// Work is split into contiguous chunks of `grain` indices that are claimed dynamically by up to
// `n_threads` threads: the calling thread plus helpers from a persistent process-wide thread pool
// (RcppParallel uses TBB's pool; spawning std::threads per call cost ~200 us per region). The body
// receives [begin, end) ranges. Results are deterministic as long as each index writes only to
// its own output slot (as in all upstream Workers), whatever the number of threads.
//
// The pool serves several regions at once (concurrent calls from different Python threads and
// nested parallel_for calls share it). The caller always works on its own region and never
// waits for helpers that have not started, so regions complete even if no helper is available
// (pool busy, thread creation failed, or after fork()). Exceptions thrown by the body are
// captured and the first one is rethrown in the caller.

#ifndef DADA2_PARALLEL_H
#define DADA2_PARALLEL_H

#include <algorithm>
#include <atomic>
#include <cstddef>
#include <exception>
#include <mutex>
#include <type_traits>

namespace dada2 {
namespace detail {

// Runs `run(ctx)` on the calling thread and on up to `n_helpers` pool threads; returns when the
// caller's run has returned and every helper that started has finished. `run` must not throw.
void run_on_pool(void (*run)(void *), void *ctx, std::size_t n_helpers);

// Upper bound on useful threads (hardware concurrency, at least 1).
std::size_t max_threads();

}  // namespace detail

template <typename F>
void parallel_for(std::size_t begin, std::size_t end, F &&body, std::size_t grain, int n_threads) {
  if (end <= begin) return;
  if (grain < 1) grain = 1;
  const std::size_t n = end - begin;
  const std::size_t nchunks = (n + grain - 1) / grain;
  if (n_threads <= 1 || nchunks <= 1) {
    body(begin, end);
    return;
  }
  // Threads beyond the number of chunks or of hardware threads cannot help.
  std::size_t nthreads = std::min<std::size_t>(static_cast<std::size_t>(n_threads), nchunks);
  nthreads = std::min(nthreads, detail::max_threads());

  struct Ctx {
    std::size_t begin, end, grain, nchunks;
    std::remove_reference_t<F> *body;
    std::atomic<std::size_t> next{0};
    std::atomic<bool> failed{false};
    std::exception_ptr error;
    std::mutex error_mutex;
  } ctx;
  ctx.begin = begin;
  ctx.end = end;
  ctx.grain = grain;
  ctx.nchunks = nchunks;
  ctx.body = &body;

  auto run = [](void *p) {
    Ctx &c = *static_cast<Ctx *>(p);
    for (;;) {
      if (c.failed.load(std::memory_order_relaxed)) return;
      const std::size_t k = c.next.fetch_add(1);
      if (k >= c.nchunks) return;
      const std::size_t b = c.begin + k * c.grain;
      const std::size_t e = std::min(c.end, b + c.grain);
      try {
        (*c.body)(b, e);
      } catch (...) {
        std::lock_guard<std::mutex> lock(c.error_mutex);
        if (!c.error) c.error = std::current_exception();
        c.failed.store(true);
        return;
      }
    }
  };
  detail::run_on_pool(run, &ctx, nthreads - 1);
  if (ctx.error) std::rethrow_exception(ctx.error);
}

}  // namespace dada2

#endif  // DADA2_PARALLEL_H
