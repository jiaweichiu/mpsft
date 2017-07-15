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
  constexpr Real snr = 10; // Rather high SNR.
  constexpr Real threshold = 0.2;

  ModeMap mm = {
    {103, Cplex(3.5, 1.1)},
    {660, Cplex(-2.4, 1.5)},
  };
  CplexArray xh = GenerateXhat(n, mm, snr);

  FFTPlan plan(n, FFTW_BACKWARD);
  CplexArray x(n);
  plan.Run(xh, &x);

  ModeMap ans_mm;
  mpsft(x, bins, delta, trials, threshold, &ans_mm);
}

} // namespace mps