#include "catch.hpp"

#include "base.h"
#include "window.h"

namespace mps {

constexpr Int n = 1109;
constexpr Int bins = 5;
constexpr Real delta = 1e-6;

TEST_CASE("WindowBasic", "") {
  Window window(n, bins, delta);

  const Int p = window.p();
  const Int p2 = window.p2();

  CplexArray a(n);
  a.Clear();
  a[0] = window.wt(0);
  for (Int i = 1; i <= p2; ++i) {
    a[i] = a[n - i] = window.wt(i);
  }

  CplexArray ah(n);
  FFTPlan plan(n, -1);
  plan.Run(a, &ah);

  for (Int i = 0; i < n; ++i) {
    REQUIRE(IM(ah[i]) == Approx(0));
  }

  REQUIRE(RE(ah[0]) == Approx(1));
  REQUIRE(RE(ah[1]) == Approx(1));
  REQUIRE(RE(ah[n - 1]) == Approx(1));
  REQUIRE(RE(ah[2]) == Approx(1));
  REQUIRE(RE(ah[n - 2]) == Approx(1));

  for (Int i = 1; i < n / (4 * bins); ++i) {
    REQUIRE(RE(ah[i]) > 0.5);
    REQUIRE(RE(ah[n - i]) > 0.5);
  }

  const Int w = n / (2 * bins);
  for (Int i = w; i < n - w; ++i) {
    REQUIRE(std::abs(RE(ah[i])) < delta);
  }

  // Check against SampleInFreq.
  REQUIRE(window.SampleInFreq(0) == Approx(RE(ah[0])));
  const Int n2 = (n - 1) / 2;
  for (Int i = 1; i <= n2; ++i) {
    const Real xi = Real(i) / n;
    REQUIRE(window.SampleInFreq(xi) == Approx(RE(ah[i])));
    REQUIRE(window.SampleInFreq(-xi) == Approx(RE(ah[n - i])));
  }
}

} // namespace mps