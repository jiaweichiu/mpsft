#include "catch.hpp"

#include "base.h"

namespace mps {

namespace {

void Fill(CplexArray *v) {
  // Just to see if it will crash.
  for (Int i = 0; i < v->Size(); ++i) {
    (*v)[i] = Cplex(0.5, -0.5);
  }
}

} // namespace

TEST_CASE("CplexArrayBasic", "") {
  CplexArray u;
  REQUIRE(u.Size() == 0);
  Fill(&u);

  u.Resize(50);
  REQUIRE(u.Size() == 50);
  Fill(&u);

  CplexArray v(100);
  REQUIRE(v.Size() == 100);
  Fill(&v);

  v.Resize(10);
  REQUIRE(v.Size() == 10);
  Fill(&v);

  v.Resize(10000);
  REQUIRE(v.Size() == 10000);
  Fill(&v);
}

TEST_CASE("FFTPlanBasic", "") {
  constexpr int n = 7;
  FFTPlan plan(n, FFTW_FORWARD, true);
  CplexArray a(n);
  Fill(&a);
  plan.RunInPlace(&a);
  REQUIRE(RE(a[0]) == Approx(3.5));
  REQUIRE(IM(a[0]) == Approx(-3.5));

  FFTPlan plan2(n, FFTW_BACKWARD, true);
  plan2.RunInPlace(&a);
  for (Int i = 0; i < a.Size(); ++i) {
    // Different from 0.5 because FFTW does unnormalized FFTs.
    REQUIRE(RE(a[i]) == Approx(3.5));
    REQUIRE(IM(a[i]) == Approx(-3.5));
  }
}

} // namespace mps
