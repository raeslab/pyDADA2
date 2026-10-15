// pydada2 runtime support replacing R/Rcpp services in the ported DADA2 core:
//   Rcpp::stop               -> dada2_stop (throws dada2::Error / std::invalid_argument)
//   Rprintf                  -> dada2_log  (optional handler of the calling thread)
//   Rcpp::checkUserInterrupt -> dada2_check_interrupt (optional handler of the calling thread)
//
// The handlers are thread-local: each thread that calls into the core installs its own (via
// ScopedHandlers), so concurrent calls from several threads do not interfere. They are only
// invoked on the thread that installed them, never inside parallel_for worker bodies.

#ifndef DADA2_RUNTIME_H
#define DADA2_RUNTIME_H

#include <atomic>
#include <functional>
#include <stdexcept>
#include <string>

namespace dada2 {

// Internal error raised by the ported core (maps to RuntimeError in Python).
class Error : public std::runtime_error {
 public:
  using std::runtime_error::runtime_error;
};

// Raised when the interrupt handler reports a pending interrupt.
class Interrupted : public std::runtime_error {
 public:
  Interrupted() : std::runtime_error("Interrupted") {}
};

using LogHandler = std::function<void(const std::string &)>;
using InterruptHandler = std::function<bool()>;

// Install handlers for the calling thread (nullptr/empty to disable).
void set_log_handler(LogHandler handler);
void set_interrupt_handler(InterruptHandler handler);

// RAII installer used by the bindings for the duration of one call: installs the handlers for the
// calling thread and restores the previous ones (of the same thread) on destruction.
class ScopedHandlers {
 public:
  ScopedHandlers(LogHandler log, InterruptHandler interrupt);
  ~ScopedHandlers();
  ScopedHandlers(const ScopedHandlers &) = delete;
  ScopedHandlers &operator=(const ScopedHandlers &) = delete;

 private:
  LogHandler prev_log_;
  InterruptHandler prev_interrupt_;
};

// RAII: makes dada2_check_interrupt() on the calling thread also throw Interrupted once *flag is
// set (used to stop concurrent work, e.g. the other samples of a batch, after an error or an
// interrupt). Restores the previous flag of the thread on destruction.
class ScopedCancel {
 public:
  explicit ScopedCancel(const std::atomic<bool> *flag);
  ~ScopedCancel();
  ScopedCancel(const ScopedCancel &) = delete;
  ScopedCancel &operator=(const ScopedCancel &) = delete;

 private:
  const std::atomic<bool> *prev_;
};

}  // namespace dada2

// printf-style helpers used by the ported C-style code.
[[noreturn]] void dada2_stop(const char *fmt, ...);
[[noreturn]] void dada2_stop_invalid(const char *fmt, ...);
void dada2_log(const char *fmt, ...);
void dada2_check_interrupt();

#endif  // DADA2_RUNTIME_H
