/*
 * Copyright (c) 2017 Jiawei Chiu
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301 USA
 *
 */
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

TEST_CASE("IntArrayBasic", "") {
  IntArray u(55);
  REQUIRE(u.size() == 55);
  for (Int i = 0; i < 55; ++i) {
    u[i] = 333;
  }
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
  REQUIRE(PosMod(tf.a_inv * tf.a, n) == 1);
}

TEST_CASE("GenerateXhatBasic", "") {
  constexpr Int n = 5; // Prime.
  constexpr double sigma = 8.0;
  constexpr Cplex coef = Cplex(0.5, 0.6);
  const ModeMap mm = {{2, coef}};

  CplexArray xh(n);
  GenerateXhat(n, mm, sigma, &xh);
  double noise_energy = 0;
  noise_energy += AbsSq(RE(xh[0]), IM(xh[0]));
  noise_energy += AbsSq(RE(xh[1]), IM(xh[1]));
  noise_energy += AbsSq(RE(xh[3]), IM(xh[3]));
  noise_energy += AbsSq(RE(xh[4]), IM(xh[4]));
  REQUIRE(std::sqrt(noise_energy) == Approx(sigma));
}

TEST_CASE("GenerateXhatFFT", "") {
  constexpr Int n = 10000; // Prime.
  constexpr double sigma = 5.5;
  const ModeMap mm;

  CplexArray xh(n);
  GenerateXhat(n, mm, sigma, &xh);
  CplexArray x(n);
  FFTPlan plan(n, FFTW_BACKWARD);
  plan.Run(xh, &x);

  // Mean square error in time domain.
  const double mse = std::sqrt(x.Energy() / x.size());
  REQUIRE(mse == Approx(sigma));
}

TEST_CASE("SinTwoPiBasic", "") {
  const Int n = 1234;
  for (Int i = 0; i <= n; ++i) {
    const double x = double(i) / double(n);
    REQUIRE(SinTwoPi(x) == Approx(std::sin(2 * M_PI * x)));
  }
}

/*TEST_CASE("SinCosTwoPiBasic", "") {
  const Int n = 1234;
  for (Int i = -n; i <= n + n; ++i) {
    const double x = double(i) / double(n);
    const Cplex z = Sinusoid(x);
    REQUIRE(IM(z) == Approx(std::sin(2 * M_PI * x)));
    REQUIRE(RE(z) == Approx(std::cos(2 * M_PI * x)));
  }
}*/

TEST_CASE("DivideBasic", "") {
  // idx selects the prime.
  for (int idx = 10; idx < 25; ++idx) {
    const Int d = kPrimes[idx];
    const Int multiplier = kPrimesMagic[idx].multiplier;
    const int shift = kPrimesMagic[idx].shift;
    for (int i = 0; i < 100; ++i) {
      const Int n = RandomInt();
      const Int ans = ApplyMagic(n, multiplier, shift);
      REQUIRE(ans == (n / d));
    }
  }
}

} // namespace mps
