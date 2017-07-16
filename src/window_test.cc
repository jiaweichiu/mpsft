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
#include "window.h"

namespace mps {

constexpr Int n = 1109;
constexpr Int bins = 5;
constexpr double delta = 1e-6;

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
    const double xi = double(i) / n;
    REQUIRE(window.SampleInFreq(xi) == Approx(RE(ah[i])));
    REQUIRE(window.SampleInFreq(-xi) == Approx(RE(ah[n - i])));
  }
}

} // namespace mps
