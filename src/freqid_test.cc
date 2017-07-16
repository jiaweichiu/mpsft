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
#include "freqid.h"

namespace mps {

TEST_CASE("SVDBasic", "") {
  Cplex u0 = Cplex(-2.3, 5.2);
  Cplex u1 = Cplex(-5.5, -3.1);
  Cplex u2 = Cplex(4.0, 3.0);

  double sigma[2];
  bool r = MatPencil(u0, u1, u2, sigma);
  REQUIRE(sigma[0] == Approx(11.34419934));
  REQUIRE(sigma[1] == Approx(0.91057194));
  REQUIRE_FALSE(r);
}

TEST_CASE("SVDBasic2", "") {
  Cplex u0 = Cplex(2.3, 5.2);
  Cplex u1 = Cplex(5.5, -3.1);
  Cplex u2 = Cplex(4.0, 3.0);

  double sigma[2];
  bool r = MatPencil(u0, u1, u2, sigma);
  REQUIRE(sigma[0] == Approx(9.80673965));
  REQUIRE(sigma[1] == Approx(5.77476038));
  REQUIRE(r);
}

TEST_CASE("SVDBasic3", "") {
  Cplex u0 = Cplex(-15.3, -5.2);
  Cplex u1 = Cplex(8.0, 3.0);
  Cplex u2 = Cplex(-1.5, 3.5);

  double sigma[2];
  bool r = MatPencil(u0, u1, u2, sigma);
  REQUIRE(sigma[0] == Approx(21.27794592));
  REQUIRE(sigma[1] == Approx(12.53032391));
  REQUIRE_FALSE(r);
}

} // namespace mps