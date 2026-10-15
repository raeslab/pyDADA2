// Statistical helpers (see stats.h). Clean-room implementation from the literature plus
// Boost.Math (BSL-1.0); no R code is used.
//
// ppois_upper(n, lambda) = P(X > n) = P(X >= m), m = n + 1, X ~ Poisson(lambda).
//
//  * lambda >= m: the result is >= ~0.45 and Boost's regularized incomplete gamma
//    gamma_p(m, lambda) is accurate to a few ulps there.
//  * lambda < m (the upper tail is the small side): Boost's gamma_p loses accuracy for large m
//    (rel. 2e-10 at m ~ 1e6) and returns exactly 0 where the tail is subnormal. Instead:
//
//        P(X >= m) = dpois(m; lambda) * S,   S = sum_{j>=0} prod_{i=1..j} lambda / (m + i)
//
//    (the textbook series of the regularized incomplete gamma / Kummer's 1F1(1; m+1; lambda)),
//    with dpois in the saddle-point form of C. Loader, "Fast and accurate computation of
//    binomial probabilities" (2000):
//
//        dpois(m; lambda) = exp(-stirlerr(m) - bd0(m, lambda)) / sqrt(2 pi m)
//        stirlerr(x)      = log(x!) - [(x + 1/2) log x - x + log sqrt(2 pi)]
//        bd0(x, mu)       = x log(x / mu) + mu - x
//
//    bd0 is evaluated with the series in v = (x - mu)/(x + mu),
//        bd0 = (x - mu) v + 2 x sum_{j>=1} v^(2j+1) / (2j + 1),
//    which has no cancellation when x ~ mu (the naive formula loses ~1e-11 relative at
//    m ~ 1e5). Results below DBL_MIN are computed as exp(log P) with a single rounding, so
//    subnormal tails are kept (R returns them as well).
//
// Measured against a 45-digit reference (mpmath) on n in [0, 2e6], lambda/n in [0.05, 2],
// lambda down to 1e-320 and a dense grid around the underflow boundary: max relative error
// 6e-13 (largest near p ~ 1e-300, where the exponent is ~700), identical zero/non-zero pattern
// as R. See tests/parity/test_core_lowlevel.py::test_ppois_upper_extended.

#include "stats.h"

#include <cmath>
#include <limits>

#include <boost/math/policies/policy.hpp>
#include <boost/math/special_functions/gamma.hpp>

namespace {

// Never throw; never promote double to long double (long double differs between platforms,
// which would make results platform dependent).
using Policy = boost::math::policies::policy<
    boost::math::policies::domain_error<boost::math::policies::errno_on_error>,
    boost::math::policies::pole_error<boost::math::policies::errno_on_error>,
    boost::math::policies::overflow_error<boost::math::policies::errno_on_error>,
    boost::math::policies::underflow_error<boost::math::policies::ignore_error>,
    boost::math::policies::evaluation_error<boost::math::policies::errno_on_error>,
    boost::math::policies::promote_double<false>>;

constexpr double kLnSqrt2Pi = 0.918938533204672741780329736406;  // log(sqrt(2 pi))
constexpr double kTwoPi = 6.283185307179586476925286766559;
constexpr double kEps = std::numeric_limits<double>::epsilon();
constexpr double kLogDblMin = -708.39641853226410622;  // log(DBL_MIN)

// Remainder of Stirling's approximation, log(x!) - [(x + 1/2) log x - x + log sqrt(2 pi)],
// for integer-valued x >= 1.
double stirlerr(double x) {
  if (x <= 15.0) {
    // x! is exact in double for x <= 22.
    double fact = 1.0;
    for (double k = 2.0; k <= x; k += 1.0) fact *= k;
    return std::log(fact) - (x + 0.5) * std::log(x) + x - kLnSqrt2Pi;
  }
  // Stirling series 1/(12x) - 1/(360x^3) + 1/(1260x^5) - 1/(1680x^7) + 1/(1188x^9)
  // - 691/(360360x^11); the first omitted term is < 1e-19 for x > 15.
  const double x2 = 1.0 / (x * x);
  double s = -691.0 / 360360.0;
  s = s * x2 + 1.0 / 1188.0;
  s = s * x2 - 1.0 / 1680.0;
  s = s * x2 + 1.0 / 1260.0;
  s = s * x2 - 1.0 / 360.0;
  s = s * x2 + 1.0 / 12.0;
  return s / x;
}

// Deviance term x log(x / mu) + mu - x (>= 0) for x > 0, mu > 0.
double bd0(double x, double mu) {
  const double d = x - mu;
  if (std::fabs(d) < 0.7 * (x + mu)) {
    const double v = d / (x + mu);
    const double v2 = v * v;
    double s = d * v;
    double ej = 2.0 * x * v;
    for (int j = 1; j < 2000; ++j) {
      ej *= v2;
      const double s1 = s + ej / (2 * j + 1);
      if (s1 == s) break;
      s = s1;
    }
    return s;
  }
  const double q = x / mu;
  const double l = std::isfinite(q) ? std::log(q) : std::log(x) - std::log(mu);
  return x * l + mu - x;
}

// S = sum_{j>=0} prod_{i=1..j} lambda / (m + i) for lambda < m (terms decrease; compensated sum).
double upper_series(double m, double lambda) {
  double s = 1.0, c = 0.0, t = 1.0;
  for (double j = 1.0; j < 1e9; j += 1.0) {
    t *= lambda / (m + j);
    const double y = s + t;
    c += (s - y) + t;
    s = y;
    if (t <= 0.5 * kEps * s) break;
  }
  return s + c;
}

}  // namespace

double ppois_upper(int n, double lambda) {
  if (std::isnan(lambda) || lambda < 0) return std::numeric_limits<double>::quiet_NaN();
  if (n < 0) return 1.0;
  if (lambda == 0) return 0.0;
  if (std::isinf(lambda)) return 1.0;
  const double m = static_cast<double>(n) + 1.0;
  // P(X > n) = P(X >= m) = regularized lower incomplete gamma P(m, lambda).
  if (lambda >= m) return boost::math::gamma_p(m, lambda, Policy());
  const double lx = -stirlerr(m) - bd0(m, lambda);       // log(dpois(m; lambda) * sqrt(2 pi m))
  const double f = upper_series(m, lambda) / std::sqrt(kTwoPi * m);
  const double lp = lx + std::log(f);                    // log P(X >= m); f <= 1, so lx >= lp
  if (lp < kLogDblMin) return std::exp(lp);              // subnormal or zero: one rounding
  return std::exp(lx) * f;
}
