// pydada2 runtime support (see runtime.h).

#include "runtime.h"

#include <cstdarg>
#include <cstdio>
#include <vector>

namespace {

dada2::LogHandler g_log;
dada2::InterruptHandler g_interrupt;

std::string vformat(const char *fmt, va_list ap) {
  va_list ap2;
  va_copy(ap2, ap);
  int n = std::vsnprintf(nullptr, 0, fmt, ap2);
  va_end(ap2);
  if (n < 0) return std::string(fmt);
  std::vector<char> buf(static_cast<size_t>(n) + 1);
  std::vsnprintf(buf.data(), buf.size(), fmt, ap);
  return std::string(buf.data(), static_cast<size_t>(n));
}

}  // namespace

namespace dada2 {

void set_log_handler(LogHandler handler) { g_log = std::move(handler); }
void set_interrupt_handler(InterruptHandler handler) { g_interrupt = std::move(handler); }

ScopedHandlers::ScopedHandlers(LogHandler log, InterruptHandler interrupt) {
  g_log = std::move(log);
  g_interrupt = std::move(interrupt);
}

ScopedHandlers::~ScopedHandlers() {
  g_log = nullptr;
  g_interrupt = nullptr;
}

}  // namespace dada2

void dada2_stop(const char *fmt, ...) {
  va_list ap;
  va_start(ap, fmt);
  std::string msg = vformat(fmt, ap);
  va_end(ap);
  throw dada2::Error(msg);
}

void dada2_stop_invalid(const char *fmt, ...) {
  va_list ap;
  va_start(ap, fmt);
  std::string msg = vformat(fmt, ap);
  va_end(ap);
  throw std::invalid_argument(msg);
}

void dada2_log(const char *fmt, ...) {
  if (!g_log) return;
  va_list ap;
  va_start(ap, fmt);
  std::string msg = vformat(fmt, ap);
  va_end(ap);
  g_log(msg);
}

void dada2_check_interrupt() {
  if (g_interrupt && g_interrupt()) throw dada2::Interrupted();
}
