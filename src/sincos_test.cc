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

#include <boost/simd/function/aligned_store.hpp>
#include <boost/simd/function/sincos.hpp>
#include <boost/simd/pack.hpp>

#include "base.h"
#include "sincos.h"

namespace bs = boost::simd;

namespace mps {

TEST_CASE("SinTwoPiBasic", "") {
  const int32_t n = 1234;
  for (int32_t i = 0; i <= n; ++i) {
    const double x = double(i) / double(n);
    REQUIRE(SinTwoPiApprox(x) == Approx(std::sin(2 * M_PI * x)));
  }
}

TEST_CASE("CosTwoPiBasic", "") {
  const int32_t n = 1234;
  for (int32_t i = 0; i <= n; ++i) {
    const double x = double(i) / double(n);
    REQUIRE(CosTwoPiApprox(x) == Approx(std::cos(2 * M_PI * x)));
  }
}

TEST_CASE("SinCosTwoPiBasic_BS", "") {
  const int32_t n = 10000;
  DoubleArray a(n);
  DoubleArray out1(n);
  DoubleArray out2(n);
  for (int i = 0; i < n; ++i) {
    a[i] = RandomNormal();
  }
  using pack_t = bs::pack<double>;
  const size_t pack_card = bs::cardinal_of<pack_t>();
  for (int32_t i = 0; i < n; i += pack_card) {
    pack_t x(bs::aligned_load<pack_t>(a.data() + i));
    auto res = bs::sincos(x);
    bs::aligned_store(res.first, out1.data() + i);
    bs::aligned_store(res.second, out2.data() + i);
  }
  for (int32_t i = 0; i < n; ++i) {
    REQUIRE(::sin(a[i]) == Approx(out1[i]));
    REQUIRE(::cos(a[i]) == Approx(out2[i]));
  }
}

} // namespace mps