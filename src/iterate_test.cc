#include "catch.hpp"

#include "base.h"
#include "binner.h"
#include "iterate.h"

namespace mps {

namespace {

std::pair<Int, Int> CountGoodBad(const ModeMap &found, const ModeMap &ans) {
  Int bad = 0;
  Int good = 0;
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

  constexpr Int n = kPrimes[12];
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
  const CplexArray xh = GenerateXhat(n, mm, sigma);

  FFTPlan plan(n, FFTW_BACKWARD);
  CplexArray x(n);
  plan.Run(xh, &x);

  ModeMap ans_mm;
  Iterate(x, opt, &ans_mm);

  // for (auto& kv : ans_mm) {
  //   LOG(INFO) << kv.first << " " << kv.second;
  // }

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

  constexpr Int n = kPrimes[20];
  constexpr Int num_modes = 1000;
  constexpr double sigma = 1e-2;

  IterateOptions opt;
  opt.bins = 2 * num_modes + 1;
  opt.window_delta = 1e-6;
  opt.trials = 3;
  opt.bin_threshold = 0.01;
  opt.window_threshold = 0.1;

  // Generate a list of random coefficients, each of magnitude 1.0.
  ModeMap mm;
  for (Int i = 0; i < num_modes; ++i) {
    Cplex coef(RandomNormal(), RandomNormal());
    coef /= std::abs(coef);
    const Int loc = RandomInt() % n;
    mm[loc] += coef;
  }
  CplexArray xh = GenerateXhat(n, mm, sigma);
  FFTPlan plan(n, FFTW_BACKWARD);
  CplexArray x(n);
  plan.Run(xh, &x); // xh is the solution.

  ModeMap found_mm;
  Iterate(x, opt, &found_mm);
  auto result = CountGoodBad(found_mm, mm);
  REQUIRE(result.first >= 300);
  REQUIRE(result.second <= 100);

  opt.bins = (num_modes - 300) * 2 + 1;
  Iterate(x, opt, &found_mm);
  result = CountGoodBad(found_mm, mm);
  REQUIRE(result.first >= 600);
  REQUIRE(result.second <= 100);
}

// Assume sigma amount of noise. Use that to calibrate IterateOptions.
// Bin threshold and SV threshold should be ~sigma^2/bins.
// For simplicity, we keep the number of bins constant.
TEST_CASE("IterateFull", "") {
  RandomSeed(123537);

  constexpr Int n = kPrimes[20];
  constexpr Int num_modes = 1000;
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
  for (Int i = 0; i < num_modes; ++i) {
    Cplex coef(RandomNormal(), RandomNormal());
    coef /= std::abs(coef);
    const Int loc = RandomInt() % n;
    mm[loc] += coef;
  }
  CplexArray xh = GenerateXhat(n, mm, sigma);
  FFTPlan plan(n, FFTW_BACKWARD);
  CplexArray x(n);
  plan.Run(xh, &x); // xh is the solution.

  ModeMap found_mm;
  Iterate(x, opt, &found_mm);
  auto result = CountGoodBad(found_mm, mm);
  REQUIRE(result.first >= 250);
  REQUIRE(result.second <= 30);
  LOG(INFO) << "Modes found: " << result.first << " good " << result.second
            << " bad";

  Iterate(x, opt, &found_mm);
  result = CountGoodBad(found_mm, mm);
  REQUIRE(result.first >= 500);
  REQUIRE(result.second <= 30);
  LOG(INFO) << "Modes found: " << result.first << " good " << result.second
            << " bad";

  Iterate(x, opt, &found_mm);
  result = CountGoodBad(found_mm, mm);
  REQUIRE(result.first >= 630);
  REQUIRE(result.second <= 30);
  LOG(INFO) << "Modes found: " << result.first << " good " << result.second
            << " bad";

  Iterate(x, opt, &found_mm);
  result = CountGoodBad(found_mm, mm);
  REQUIRE(result.first >= 730);
  REQUIRE(result.second <= 30);
  LOG(INFO) << "Modes found: " << result.first << " good " << result.second
            << " bad";

  Iterate(x, opt, &found_mm);
  result = CountGoodBad(found_mm, mm);
  REQUIRE(result.first >= 810);
  REQUIRE(result.second <= 30);
  LOG(INFO) << "Modes found: " << result.first << " good " << result.second
            << " bad";
}

} // namespace mps