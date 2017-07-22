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
#include "binner.h"
#include "gen.h"
#include "iterate.h"

namespace mps {

namespace {

std::pair<int32_t, int32_t> CountGoodBad(const ModeMap &found,
                                         const ModeMap &ans) {
  int32_t bad = 0;
  int32_t good = 0;
  for (const auto &kv : found) {
    auto it = ans.find(kv.first);
    if (it == ans.end()) {
      ++bad;
    } else {
      ++good;
    }
  }
  return std::make_pair(good, bad);
}

} // namespace

TEST_CASE("IterateOnce", "") {
  RandomSeed(123537);

  constexpr int32_t n = kPrimes[12];
  constexpr double sigma = 0.1;

  IterateOptions opt;
  opt.bins = 11;
  opt.window_delta = 1e-6;
  opt.trials = 5;
  opt.bin_threshold = 0.2;
  opt.window_threshold = 0.1;

  ModeMap mm = {
      {313, Cplex(3.5, 1.1)}, {660, Cplex(-2.4, 1.5)},
  };
  CplexArray xh(n);
  GenerateXhat(n, mm, sigma, &xh);

  FFTPlan plan(n, FFTW_BACKWARD);
  CplexArray x(n);
  plan.Run(xh, &x);

  ModeMap ans_mm;
  Iterate(x, opt, &ans_mm);

  REQUIRE(ans_mm.size() == 2);
  auto it1 = ans_mm.find(313);
  REQUIRE(it1 != ans_mm.end());
  REQUIRE(std::abs(it1->second - Cplex(3.5, 1.1)) < 1e-1);

  auto it2 = ans_mm.find(660);
  REQUIRE(it2 != ans_mm.end());
  REQUIRE(std::abs(it2->second - Cplex(-2.4, 1.5)) < 1e-1);
}

// Use a larger size array. Use more modes.
TEST_CASE("IterateMore", "") {
  RandomSeed(123537);

  constexpr int32_t n = kPrimes[20];
  constexpr int32_t num_modes = 1000;
  constexpr double sigma = 1e-2;

  IterateOptions opt;
  opt.bins = 2 * num_modes + 1;
  opt.window_delta = 1e-6;
  opt.trials = 3;
  opt.bin_threshold = 0.01;
  opt.window_threshold = 0.1;

  // Generate a list of random coefficients, each of magnitude 1.0.
  ModeMap mm;
  GenerateModeMap(n, num_modes, &mm);

  CplexArray xh(n);
  GenerateXhat(n, mm, sigma, &xh);
  FFTPlan plan(n, FFTW_BACKWARD);
  CplexArray x(n);
  plan.Run(xh, &x); // xh is the solution.

  ModeMap found_mm;
  Iterate(x, opt, &found_mm);
  auto result = CountGoodBad(found_mm, mm);
  LOG(INFO) << "Good:" << result.first << " Bad:" << result.second;
  REQUIRE(result.first >= 300);
  REQUIRE(result.second <= 100);

  opt.bins = (num_modes - 300) * 2 + 1;
  Iterate(x, opt, &found_mm);
  result = CountGoodBad(found_mm, mm);
  LOG(INFO) << "Good:" << result.first << " Bad:" << result.second;
  REQUIRE(result.first >= 600);
  REQUIRE(result.second <= 100);
}

// Assume sigma amount of noise. Use that to calibrate IterateOptions.
// Bin threshold and SV threshold should be ~sigma^2/bins.
// For simplicity, we keep the number of bins constant.
TEST_CASE("IterateFull", "") {
  RandomSeed(123537);

  constexpr int32_t n = kPrimes[20];
  constexpr int32_t num_modes = 1000;
  constexpr double sigma = 1e-7;

  IterateOptions opt;
  opt.bins = 2 * num_modes + 1;
  opt.window_delta = 1e-6;
  opt.trials = 3;
  opt.bin_threshold = 3.0 * sigma / std::sqrt(double(opt.bins));
  opt.sv_threshold = 1.0 * sigma / std::sqrt(double(opt.bins));
  opt.window_threshold = 0.1;

  // Generate a list of random coefficients, each of magnitude 1.0.
  ModeMap mm;
  GenerateModeMap(n, num_modes, &mm);

  CplexArray xh(n);
  GenerateXhat(n, mm, sigma, &xh);
  FFTPlan plan(n, FFTW_BACKWARD);
  CplexArray x(n);
  plan.Run(xh, &x); // xh is the solution.

  ModeMap found_mm;
  Iterate(x, opt, &found_mm);
  auto result = CountGoodBad(found_mm, mm);
  LOG(INFO) << "Good:" << result.first << " Bad:" << result.second;
  REQUIRE(result.first >= 250);
  REQUIRE(result.second <= 30);

  Iterate(x, opt, &found_mm);
  result = CountGoodBad(found_mm, mm);
  LOG(INFO) << "Good:" << result.first << " Bad:" << result.second;
  REQUIRE(result.first >= 450);
  REQUIRE(result.second <= 30);

  Iterate(x, opt, &found_mm);
  result = CountGoodBad(found_mm, mm);
  LOG(INFO) << "Good:" << result.first << " Bad:" << result.second;
  REQUIRE(result.first >= 600);
  REQUIRE(result.second <= 30);

  Iterate(x, opt, &found_mm);
  result = CountGoodBad(found_mm, mm);
  LOG(INFO) << "Good:" << result.first << " Bad:" << result.second;
  REQUIRE(result.first >= 700);
  REQUIRE(result.second <= 30);
}

} // namespace mps