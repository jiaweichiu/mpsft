#include "catch.hpp"

#include "base.h"
#include "binner.h"
#include "window.h"

namespace mps {

TEST_CASE("BinnerBasic", "") {
  const Int n = 1109;
  const Int bins = 5;
  const Int bits = 2;
  const Int q = 106;
  Window win(n, bins, 1e-6);

  // Prepare x_hat and x.
  const ModeMap mm = {{550, Cplex(2.0, 0)}};
  const CplexArray x = EvaluateModes(n, mm);

  // Prepare binning.
  Transform tf(n, 3, 847, 45);
  Binner binner(win, tf, bits);

  // BinInTime.
  CplexMatrix out_time(1 + 2 * bits, bins);
  binner.BinInTime(x, q, &out_time);

  REQUIRE(RE(out_time[0][1]) == Approx(1.12652));
  REQUIRE(IM(out_time[0][1]) == Approx(-0.108838));
  for (Int i = 0; i < bins; ++i) {
    if (i != 1) {
      REQUIRE(std::abs((out_time[0])[i]) == Approx(0));
    }
  }

  // BinInFreq.
  CplexMatrix out_freq(1 + 2 * bits, bins);
  out_freq.Clear();
  binner.BinInFreq(mm, q, &out_freq); // Subtract.

  REQUIRE(RE(out_freq[0][1]) == Approx(-1.12652));
  REQUIRE(IM(out_freq[0][1]) == Approx(0.108838));
  for (Int i = 0; i < bins; ++i) {
    if (i != 1) {
      REQUIRE(std::abs(out_freq[0][i]) == Approx(0));
    }
  }

  // Compare out_time and out_freq.
  for (Int i = 0; i < 1 + 2 * bits; ++i) {
    for (Int j = 0; j < bins; ++j) {
      REQUIRE(std::abs(out_time[i][j] + out_freq[i][j]) == Approx(0));
    }
  }
}

TEST_CASE("BinnerBigger", "") {
  const Int n = kPrimes[20];
  const Int bins = 5;
  const Int bits = 2;
  const Int q = 0xFFFFFF;
  Window win(n, bins, 1e-6);

  // Prepare x_hat and x.
  const ModeMap mm = {{565336, Cplex(2.0, 1.0)}};
  const CplexArray x = EvaluateModes(n, mm);

  // Prepare binning.
  Transform tf(n, 0x3FFFFFFF, 0xEEEEEEEE, 0xDDDDDD);
  Binner binner(win, tf, bits);

  CplexMatrix out_time(1 + 2 * bits, bins);
  binner.BinInTime(x, q, &out_time);

  // BinInFreq.
  CplexMatrix out_freq(1 + 2 * bits, bins);
  out_freq.Clear();
  binner.BinInFreq(mm, q, &out_freq); // Subtract.

  for (Int i = 0; i < bins; ++i) {
    REQUIRE(std::abs(out_time[0][i] + out_freq[0][i]) == Approx(0));
  }
}

} // namespace mps