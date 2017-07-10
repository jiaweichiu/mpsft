#include "catch.hpp"

#include "base.h"

namespace mps {

TEST_CASE("CplexArrayBasic", "") {
  CplexArray u;
  REQUIRE(u.Size() == 0);
  u.Fill(Cplex(0.5, -0.5));

  u.Resize(50);
  REQUIRE(u.Size() == 50);
  u.Fill(Cplex(0.5, -0.5));

  CplexArray v(100);
  REQUIRE(v.Size() == 100);
  v.Fill(Cplex(0.5, -0.5));

  v.Resize(10);
  REQUIRE(v.Size() == 10);
  v.Fill(Cplex(0.5, -0.5));

  v.Resize(10000);
  REQUIRE(v.Size() == 10000);
  v.Fill(Cplex(0.5, -0.5));
}

TEST_CASE("FFTPlanBasic", "") {
  constexpr int n = 7;
  FFTPlan plan(n, FFTW_FORWARD, true);
  CplexArray a(n);
  a.Fill(Cplex(0.5, -0.5));
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

TEST_CASE("TransformBasic", "") {
  constexpr Int n = 536870909;  // Prime.
  Transform tf(n, 10000000, 10000001, 10000002);
  REQUIRE(tf.a == 10000000);
  REQUIRE(tf.b == 10000001);
  REQUIRE(tf.c == 10000002);
  REQUIRE(PosMod(Long(tf.a_inv) * Long(tf.a), n) == 1);
}

}  // namespace mps
