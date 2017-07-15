#include "catch.hpp"

#include "base.h"

namespace mps {

TEST_CASE("CplexArrayBasic", "") {
  CplexArray u;
  REQUIRE(u.size() == 0);
  u.Fill(Cplex(0.5, -0.5));

  u.Resize(50);
  REQUIRE(u.size() == 50);
  u.Fill(Cplex(0.5, -0.5));

  u.Reset();
  u.Reset();
  u.Reset();

  CplexArray v(100);
  REQUIRE(v.size() == 100);
  v.Fill(Cplex(0.5, -0.5));

  v.Resize(10);
  REQUIRE(v.size() == 10);
  v.Fill(Cplex(0.5, -0.5));

  v.Resize(10000);
  REQUIRE(v.size() == 10000);
  v.Fill(Cplex(0.5, -0.5));

  CplexArray w = {{0.5, 0.6}, {1.5, 1.6}};
  REQUIRE(std::abs(w[0] - Cplex(0.5, 0.6)) == Approx(0));
  REQUIRE(std::abs(w[1] - Cplex(1.5, 1.6)) == Approx(0));
}

TEST_CASE("FFTPlanBasic", "") {
  constexpr int n = 7;
  FFTPlan plan(n, FFTW_FORWARD);
  CplexArray a(n);
  CplexArray b(n);
  a.Fill(Cplex(0.5, -0.5));

  plan.Run(a, &b);
  REQUIRE(RE(b[0]) == Approx(3.5));
  REQUIRE(IM(b[0]) == Approx(-3.5));

  FFTPlan plan2(n, FFTW_BACKWARD);
  plan2.Run(b, &a);
  for (Int i = 0; i < a.size(); ++i) {
    // Different from 0.5 because FFTW does unnormalized FFTs.
    REQUIRE(RE(a[i]) == Approx(3.5));
    REQUIRE(IM(a[i]) == Approx(-3.5));
  }
}

TEST_CASE("TransformBasic", "") {
  constexpr Int n = 536870909; // Prime.
  Transform tf(n, 10000000, 10000001, 10000002);
  REQUIRE(tf.a == 10000000);
  REQUIRE(tf.b == 10000001);
  REQUIRE(tf.c == 10000002);
  REQUIRE(PosMod(Long(tf.a_inv) * Long(tf.a), n) == 1);
}

TEST_CASE("GenerateXhat", "") {
  constexpr Int n = 5;    // Prime.
  constexpr Real snr = 8; // Expected snr.
  constexpr Cplex coef = Cplex(0.5, 0.6);
  const ModeMap mm = {{2, coef}};
  const CplexArray xh = GenerateXhat(n, mm, snr);
  const Real signal_energy = AbsSq(coef);
  const Real noise_energy =
      AbsSq(xh[0]) + AbsSq(xh[1]) + AbsSq(xh[3]) + AbsSq(xh[4]);
  const Real snr2 = std::sqrt(signal_energy / noise_energy); // Measured snr.
  REQUIRE(snr2 == Approx(snr));
}

} // namespace mps
