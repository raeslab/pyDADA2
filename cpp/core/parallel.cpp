// Persistent thread pool behind dada2::parallel_for (see parallel.h).

#include "parallel.h"

#include <condition_variable>
#include <deque>
#include <memory>
#include <system_error>
#include <thread>
#include <vector>

#if defined(__unix__) || defined(__APPLE__)
#include <pthread.h>
#define DADA2_HAVE_PTHREAD_ATFORK 1
#endif

namespace dada2 {
namespace detail {
namespace {

// One parallel region. Helpers hold a shared_ptr, so the bookkeeping outlives the caller's wait.
struct Region {
  void (*run)(void *) = nullptr;
  void *ctx = nullptr;
  int running = 0;  // helpers executing run(); guarded by Pool::m_
  std::condition_variable done;
};

class Pool {
 public:
  void run(void (*fn)(void *), void *ctx, std::size_t n_helpers) {
    auto region = std::make_shared<Region>();
    region->run = fn;
    region->ctx = ctx;
    std::size_t posted = 0;
    if (n_helpers > 0) {
      std::lock_guard<std::mutex> lk(m_);
      grow(n_helpers);
      posted = std::min(n_helpers, workers_.size());
      try {
        for (std::size_t i = 0; i < posted; ++i) queue_.push_back(region);
      } catch (...) {  // bad_alloc: never leave entries pointing at this stack frame behind
        queue_.erase(std::remove(queue_.begin(), queue_.end(), region), queue_.end());
        posted = 0;
      }
    }
    for (std::size_t i = 0; i < posted; ++i) cv_.notify_one();

    fn(ctx);  // the caller works on its own region until no chunk is left

    if (posted == 0) return;
    std::unique_lock<std::mutex> lk(m_);
    // Helpers that have not started are no longer needed (all chunks are claimed).
    queue_.erase(std::remove(queue_.begin(), queue_.end(), region), queue_.end());
    region->done.wait(lk, [&] { return region->running == 0; });
  }

 private:
  // Grow the pool to n workers (called with m_ held). A failure to create a thread is not fatal:
  // regions then run with the threads that exist (or on the caller alone).
  void grow(std::size_t n) {
    n = std::min(n, max_threads() > 1 ? max_threads() - 1 : std::size_t(0));
    while (workers_.size() < n && !grow_failed_) {
      try {
        workers_.emplace_back([this] { loop(); });
      } catch (const std::system_error &) {
        grow_failed_ = true;
      } catch (const std::bad_alloc &) {
        grow_failed_ = true;
      }
    }
  }

  void loop() {
    std::unique_lock<std::mutex> lk(m_);
    for (;;) {
      cv_.wait(lk, [&] { return !queue_.empty(); });
      std::shared_ptr<Region> region = std::move(queue_.front());
      queue_.pop_front();
      ++region->running;
      lk.unlock();
      region->run(region->ctx);
      lk.lock();
      if (--region->running == 0) region->done.notify_all();
    }
  }

  std::mutex m_;
  std::condition_variable cv_;
  std::deque<std::shared_ptr<Region>> queue_;
  std::vector<std::thread> workers_;  // never joined: the pool lives until process exit
  bool grow_failed_ = false;
};

// The pool is intentionally leaked (worker threads block in wait() until the process exits;
// destroying joinable std::threads at exit would call std::terminate).
std::atomic<Pool *> g_pool{nullptr};
std::once_flag g_pool_once;

#ifdef DADA2_HAVE_PTHREAD_ATFORK
// After fork() only the forking thread exists in the child: start over with an empty pool (the
// parent's pool object, possibly with a locked mutex, is abandoned).
void reset_pool_in_child() { g_pool.store(new Pool()); }
#endif

Pool &pool() {
  std::call_once(g_pool_once, [] {
    g_pool.store(new Pool());
#ifdef DADA2_HAVE_PTHREAD_ATFORK
    pthread_atfork(nullptr, nullptr, reset_pool_in_child);
#endif
  });
  return *g_pool.load();
}

}  // namespace

void run_on_pool(void (*run)(void *), void *ctx, std::size_t n_helpers) {
  if (n_helpers == 0) {
    run(ctx);
    return;
  }
  pool().run(run, ctx, n_helpers);
}

std::size_t max_threads() {
  static const std::size_t n = [] {
    const unsigned hw = std::thread::hardware_concurrency();
    return hw > 0 ? static_cast<std::size_t>(hw) : std::size_t(1024);
  }();
  return n;
}

}  // namespace detail
}  // namespace dada2
