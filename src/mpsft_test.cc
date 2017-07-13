#include "catch.hpp"

#include "base.h"
#include "binner.h"
#include "mpsft.h"

namespace mps {

TEST_CASE("MPSFTBasic", "") {
  RandomSeed(123537);

  constexpr Int n = 1109;
  constexpr Int bins = 5;
  constexpr Real delta = 1e-6;
  constexpr Int trials = 5;
  constexpr Real snr = 15; // Rather high SNR.
  constexpr Real threshold = 0.2;

  CplexArray coef = {{3.5, 0}, {-2.4, 0}};
  vector<Int> loc = {103, 660};

  CplexArray xh = GenerateXhat(n, coef, loc, snr);

  FFTPlan plan(n, FFTW_BACKWARD);
  CplexArray x(n);
  plan.Run(xh, &x);

  mpsft(x, bins, delta, trials, threshold);
}

} // namespace mps