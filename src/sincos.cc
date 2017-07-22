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
#include "sincos.h"

namespace mps {

double SinTwoPiApprox(double a) {
  const double y = (a >= 0.75) ? -1.0 + a : (a >= 0.25) ? 0.5 - a : a;
  const double x = 4.0 * y;
  constexpr double coef1 = 1.1336481778117871;
  constexpr double coef3 = -0.13807177660048911;
  constexpr double coef5 = 0.0044907175846143066;
  constexpr double coef7 = -6.8290405376023045e-05;
  constexpr double z0 = 1.0;
  const double z1 = x;
  const double z2 = 2.0 * x * z1 - z0;
  const double z3 = 2.0 * x * z2 - z1;
  const double z4 = 2.0 * x * z3 - z2;
  const double z5 = 2.0 * x * z4 - z3;
  const double z6 = 2.0 * x * z5 - z4;
  const double z7 = 2.0 * x * z6 - z5;
  return coef1 * z1 + coef3 * z3 + coef5 * z5 + coef7 * z7;
}

// cos(x) = sin(x+pi/2).
// cos(2*pi*x) = sin(2*pi(x+0.25)).
// 4th quadrant: x -> x+0.25-1.0 = x-0.75.
// 3rd quadrant: x -> -1.0+(x+0.25) = x-0.75
// 2nd quadrant: x -> 0.5-(x+0.25) = 0.25-x.
// 1st quadrant: x -> 0.5-(x+0.25) = 0.25-x.
double CosTwoPiApprox(double a) {
  const double y = (a >= 0.50) ? a - 0.75 : 0.25 - a;
  const double x = 4.0 * y;
  constexpr double coef1 = 1.1336481778117871;
  constexpr double coef3 = -0.13807177660048911;
  constexpr double coef5 = 0.0044907175846143066;
  constexpr double coef7 = -6.8290405376023045e-05;
  constexpr double z0 = 1.0;
  const double z1 = x;
  const double z2 = 2.0 * x * z1 - z0;
  const double z3 = 2.0 * x * z2 - z1;
  const double z4 = 2.0 * x * z3 - z2;
  const double z5 = 2.0 * x * z4 - z3;
  const double z6 = 2.0 * x * z5 - z4;
  const double z7 = 2.0 * x * z6 - z5;
  return coef1 * z1 + coef3 * z3 + coef5 * z5 + coef7 * z7;
}

} // namespace mps
