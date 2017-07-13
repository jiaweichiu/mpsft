#include "catch.hpp"

#include "base.h"
#include "binner.h"
#include "window.h"

namespace mps {

TEST_CASE("BinnerBasic", "") {
  const Int n = 1109;
  const Int bins = 5;
  Window win(n, bins, 1e-6);

  // Prepare x_hat and x.
  const CplexArray coef = {{2.0, 0}};
  const vector<Int> loc = {550};
  const CplexArray x = EvaluateModes(n, coef, loc);

  // Prepare binning.
  Transform tf(n, 3, 847, 45);
  TauSet taus;
  taus.q = 106;
  taus.list_s = {3, 7};

  // BinInTime.
  CplexArray scratch(bins);
  CplexMatrix out_time(taus.size(), bins);

  FFTPlan plan(bins, -1);
  BinInTime(win, tf, taus, x, &plan, &out_time, &scratch);
  REQUIRE(RE(out_time[0][1]) == Approx(1.12652));
  REQUIRE(IM(out_time[0][1]) == Approx(-0.108838));
  for (Int i = 0; i < bins; ++i) {
    if (i != 1) {
      REQUIRE(std::abs((out_time[0])[i]) == Approx(0));
    }
  }

  // BinInFreq.
  CplexMatrix out_freq(taus.size(), bins);
  out_freq.Clear();
  BinInFreq(win, tf, taus, coef, loc, &out_freq); // Subtract.

  REQUIRE(RE(out_freq[0][1]) == Approx(-1.12652));
  REQUIRE(IM(out_freq[0][1]) == Approx(0.108838));
  for (Int i = 0; i < bins; ++i) {
    if (i != 1) {
      REQUIRE(std::abs(out_freq[0][i]) == Approx(0));
    }
  }

  // Compare out_time and out_freq.
  for (Int i = 0; i < taus.size(); ++i) {
    for (Int j = 0; j < bins; ++j) {
      REQUIRE(std::abs(out_time[i][j] + out_freq[i][j]) == Approx(0));
    }
  }
}

} // namespace mps