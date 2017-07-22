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
#include "integer.h"

namespace mps {

TEST_CASE("MulMod", "") {
  const int64_t divisor = kPrimes[15];
  for (int i = 0; i < 1000; ++i) {
    const int64_t a = RandomInt32();
    const int64_t b = RandomInt32();
    REQUIRE(MulMod(a, b, divisor) == (a * b) % divisor);
  }
}

} // namespace mps