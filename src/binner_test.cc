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
  const ModeMap mm = {{550, Cplex(2.0, 0)}};
  const CplexArray x = EvaluateModes(n, mm);

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
  BinInFreq(win, tf, taus, mm, &out_freq); // Subtract.

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

TEST_CASE("BinnerBigger", "") {
  const Int n = kPrimes[20];
  const Int bins = 5;
  Window win(n, bins, 1e-6);

  // Prepare x_hat and x.
  const ModeMap mm = {{565336, Cplex(2.0, 1.0)}};
  const CplexArray x = EvaluateModes(n, mm);

  // Prepare binning.
  Transform tf(n, 0x3FFFFFFF, 0xEEEEEEEE, 0xDDDDDD);
  TauSet taus;
  taus.q = 0xFFFFFF;

  CplexArray scratch(bins);
  CplexMatrix out_time(taus.size(), bins);
  FFTPlan plan(bins, -1);
  BinInTime(win, tf, taus, x, &plan, &out_time, &scratch);

  for (Int i = 0; i < bins; ++i) {
    LOG(INFO) << out_time[0][i];
  }

  CplexMatrix out_freq(taus.size(), bins);
  out_freq.Clear();
  BinInFreq(win, tf, taus, mm, &out_freq); // Subtract.
  for (Int i = 0; i < bins; ++i) {
    LOG(INFO) << out_freq[0][i];
  }

  for (Int i = 0; i < bins; ++i) {
    REQUIRE(std::abs(out_time[0][i] + out_freq[0][i]) == Approx(0));
  }
}

} // namespace mps