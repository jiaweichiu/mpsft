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

TEST_CASE("CplexArray", "") {
  CplexArray u;
  REQUIRE(u.size() == 0);
  u.fill(Cplex(0.5, -0.5));

  u.resize(50);
  REQUIRE(u.size() == 50);
  u.fill(Cplex(0.5, -0.5));

  u.reset();
  u.reset();
  u.reset();

  CplexArray v(100);
  REQUIRE(v.size() == 100);
  v.fill(Cplex(0.5, -0.5));

  v.resize(10);
  REQUIRE(v.size() == 10);
  v.fill(Cplex(0.5, -0.5));

  v.resize(10000);
  REQUIRE(v.size() == 10000);
  v.fill(Cplex(0.5, -0.5));

  CplexArray w = {{0.5, 0.6}, {1.5, 1.6}};
  REQUIRE(std::abs(w[0] - Cplex(0.5, 0.6)) == Approx(0));
  REQUIRE(std::abs(w[1] - Cplex(1.5, 1.6)) == Approx(0));
}

TEST_CASE("IntArray", "") {
  Int32Array u(55);
  REQUIRE(u.size() == 55);
  for (int32_t i = 0; i < 55; ++i) {
    u[i] = 333;
  }
  for (int32_t i = 0; i < 55; ++i) {
    REQUIRE(u[i] == 333);
  }
}

TEST_CASE("FFTPlan", "") {
  constexpr int n = 7;
  FFTPlan plan(n, FFTW_FORWARD);
  CplexArray a(n);
  CplexArray b(n);
  a.fill(Cplex(0.5, -0.5));

  plan.Run(a, &b);
  REQUIRE(RE(b[0]) == Approx(3.5));
  REQUIRE(IM(b[0]) == Approx(-3.5));

  FFTPlan plan2(n, FFTW_BACKWARD);
  plan2.Run(b, &a);
  for (int32_t i = 0; i < a.size(); ++i) {
    // Different from 0.5 because FFTW does unnormalized FFTs.
    REQUIRE(RE(a[i]) == Approx(3.5));
    REQUIRE(IM(a[i]) == Approx(-3.5));
  }
}

TEST_CASE("DivideMagic", "") {
  // idx selects the prime.
  for (int idx = 10; idx < 25; ++idx) {
    const int32_t d = kPrimes[idx];
    const int32_t multiplier = kPrimesMagic[idx].multiplier;
    const int shift = kPrimesMagic[idx].shift;
    for (int i = 0; i < 100; ++i) {
      const int32_t n = RandomInt32();
      const int32_t ans = ApplyMagic(n, multiplier, shift);
      REQUIRE(ans == (n / d));
    }
  }
}

} // namespace mps
