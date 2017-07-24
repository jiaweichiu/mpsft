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
#include "gen.h"

namespace mps {

TEST_CASE("GenerateXhatBasic", "") {
  constexpr int32_t n = 5; // Prime.
  constexpr double sigma = 8.0;
  constexpr Cplex coef(0.5, 0.6);
  const ModeMap mm = {{2, coef}};

  CplexArray xh(n);
  GenerateXhat(n, mm, sigma, &xh);
  xh[2] = 0;
  REQUIRE(std::sqrt(xh.energy()) == Approx(sigma));
}

TEST_CASE("GenerateXhatFFT", "") {
  constexpr int32_t n = 10000; // Prime.
  constexpr double sigma = 5.5;
  const ModeMap mm;

  CplexArray xh(n);
  GenerateXhat(n, mm, sigma, &xh);
  CplexArray x(n);
  FFTPlan plan(n, FFTW_BACKWARD);
  plan.Run(xh, &x);

  // Mean square error in time domain.
  const double mse = std::sqrt(x.energy() / x.size());
  REQUIRE(mse == Approx(sigma));
}

TEST_CASE("GenerateXhatAltBasic", "") {
  constexpr int32_t n = 5; // Prime.
  constexpr double sigma = 8.0;
  constexpr Cplex coef(0.5, 0.6);
  const ModeMap mm = {{2, coef}};

  double total_energy = 0;
  CplexArray xh(n);
  constexpr int num_trials = 100000;
  for (int i = 0; i < num_trials; ++i) {
    GenerateXhatAlt(n, mm, sigma, &xh);
    xh[2] -= coef;
    total_energy += xh.energy();
  }
  const double noise_energy = total_energy / num_trials;
  REQUIRE(std::abs(noise_energy - sigma * sigma) < 0.1);
}

} // namespace mps
