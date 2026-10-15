// pydada2 runtime support replacing R/Rcpp services in the ported DADA2 core:
//   Rcpp::stop               -> dada2_stop (throws dada2::Error / std::invalid_argument)
//   Rprintf                  -> dada2_log  (optional handler, main thread only)
//   Rcpp::checkUserInterrupt -> dada2_check_interrupt (optional handler, main thread only)

#ifndef DADA2_RUNTIME_H
#define DADA2_RUNTIME_H

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

// Install handlers (nullptr/empty to disable). Not thread-safe; set from the main thread.
void set_log_handler(LogHandler handler);
void set_interrupt_handler(InterruptHandler handler);

// RAII installer used by the bindings for the duration of one call.
class ScopedHandlers {
 public:
  ScopedHandlers(LogHandler log, InterruptHandler interrupt);
  ~ScopedHandlers();
  ScopedHandlers(const ScopedHandlers &) = delete;
  ScopedHandlers &operator=(const ScopedHandlers &) = delete;
};

}  // namespace dada2

// printf-style helpers used by the ported C-style code.
[[noreturn]] void dada2_stop(const char *fmt, ...);
[[noreturn]] void dada2_stop_invalid(const char *fmt, ...);
void dada2_log(const char *fmt, ...);
void dada2_check_interrupt();

#endif  // DADA2_RUNTIME_H
