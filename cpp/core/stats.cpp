// Statistical helpers (see stats.h). Implemented with Boost.Math (BSL-1.0); no R code is used.

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

}  // namespace

double ppois_upper(int n, double lambda) {
  if (std::isnan(lambda) || lambda < 0) return std::numeric_limits<double>::quiet_NaN();
  if (n < 0) return 1.0;
  if (lambda == 0) return 0.0;
  if (std::isinf(lambda)) return 1.0;
  // P(X > n) = P(X >= n+1) = regularized lower incomplete gamma P(n+1, lambda).
  return boost::math::gamma_p(static_cast<double>(n) + 1.0, lambda, Policy());
}
