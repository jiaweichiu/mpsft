#include "catch.hpp"

#include "base.h"
#include "binner.h"
#include "iterate.h"

namespace mps {

namespace {

std::pair<Int, Int> CountGoodBad(const ModeMap& found, const ModeMap& ans) {
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

}  // namespace

/*TEST_CASE("IterateBasic", "") {
  RandomSeed(123537);

  constexpr Int n = 1109;
  constexpr double sigma = 0.2;

  IterateOptions opt;
  opt.bins = 5;
  opt.window_delta = 1e-6;
  opt.trials = 5;
  opt.bin_threshold = 0.2;
  opt.window_threshold = 0.1;

  ModeMap mm = {
      {103, Cplex(3.5, 1.1)}, {660, Cplex(-2.4, 1.5)},
  };
  const CplexArray xh = GenerateXhat(n, mm, sigma);

  FFTPlan plan(n, FFTW_BACKWARD);
  CplexArray x(n);
  plan.Run(xh, &x);

  ModeMap ans_mm;
  Iterate(x, opt, &ans_mm);

  REQUIRE(ans_mm.size() == 2);
  auto it1 = ans_mm.find(103);
  REQUIRE(it1 != ans_mm.end());
  REQUIRE(std::abs(it1->second - Cplex(3.5, 1.1)) < 1e-1);

  auto it2 = ans_mm.find(660);
  REQUIRE(it2 != ans_mm.end());
  REQUIRE(std::abs(it2->second - Cplex(-2.4, 1.5)) < 1e-1);
}*/


TEST_CASE("IterateMore", "") {
  RandomSeed(123537);

  constexpr Int n = kPrimes[20];
  constexpr Int num_modes = 1000;
  constexpr double sigma = 1e-7;

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

} // namespace mps