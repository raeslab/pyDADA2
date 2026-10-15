// Minimal Catch2 tests of the C++ core API (cpp/core/api.h, stats.h), independent of Python.
// The R-parity checks live in the Python test suite (tests/parity); these pin basic invariants so
// the core can be checked on its own (cmake -DDADA2_BUILD_CPP_TESTS=ON, then ctest).

#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

#include "api.h"
#include "stats.h"

TEST_CASE("nwalign aligns identical sequences without gaps", "[evaluate]") {
  const std::string s = "ACGTACGTTTGACCAGGT";
  auto al = dada2::nwalign(s, s, 5, -4, -8, -8, -1, true);
  REQUIRE(al.first == s);
  REQUIRE(al.second == s);
  auto ev = dada2::eval_pair(al.first, al.second);
  REQUIRE(ev[0] == static_cast<int>(s.size()));
  REQUIRE(ev[1] == 0);
  REQUIRE(ev[2] == 0);
}

TEST_CASE("nwalign places a single internal deletion", "[evaluate]") {
  auto al = dada2::nwalign("ACGTACGT", "ACGACGT", 5, -4, -8, -8, -1, true);
  REQUIRE(al.first == "ACGTACGT");
  REQUIRE(al.second == "ACG-ACGT");
  auto ev = dada2::eval_pair(al.first, al.second);
  REQUIRE(ev == std::array<int, 3>{7, 0, 1});
  REQUIRE(dada2::pair_consensus(al.first, al.second, 1, false) == "ACGTACGT");
}

TEST_CASE("eval_pair rejects alignments of different lengths", "[evaluate]") {
  REQUIRE_THROWS_AS(dada2::eval_pair("AC-", "AC"), std::invalid_argument);
}

TEST_CASE("is_acgt", "[evaluate]") {
  // As R's C_isACGT: lowercase is not ACGT, the empty string is.
  auto r = dada2::is_acgt({"ACGT", "ACGN", "acgt", ""});
  REQUIRE(r == std::vector<bool>{true, false, false, true});
}

TEST_CASE("ppois_upper edge cases", "[stats]") {
  REQUIRE(ppois_upper(-1, 3.0) == 1.0);
  REQUIRE(ppois_upper(3, 0.0) == 0.0);
  // P(X > 0) = 1 - exp(-lambda)
  REQUIRE(std::fabs(ppois_upper(0, 2.0) - (1.0 - std::exp(-2.0))) < 1e-15);
}
