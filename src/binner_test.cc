#include "catch.hpp"

#include "base.h"
#include "binner.h"
#include "window.h"

namespace mps {

TEST_CASE("BinnerBasic", "") {
  const Int n = 1109;
  const Int bins = 5;
  Window win(n, bins, 1e-6);

  // Prepare x_hat.
  CplexArray coef(1);
  coef[0] = Cplex(2.0, 0);
  vector<Int> loc(1);
  loc[0] = 550;

  // Prepare x.
  CplexArray x(n);
  x.Clear();
  for (Int t = 0; t < n; ++t) {
    x[t] +=
        coef[0] * Sinusoid((2.0 * M_PI) * Mod(loc[0] * Real(t), n) / Real(n));
  }

  // Prepare binning.
  Transform tf(n, 3, 847, 45);
  TauSet taus;
  taus.q = 106;

  // BinInTime.
  CplexArray scratch(bins);
  CplexMatrix out_time(1, bins);
  
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
  CplexArray out(bins);
  out.Clear();
  BinInFreq(win, tf, taus, coef, loc, &out); // Subtract.

  REQUIRE(RE(out[1]) == Approx(-1.12652));
  REQUIRE(IM(out[1]) == Approx(0.108838));
  for (Int i = 0; i < bins; ++i) {
    if (i != 1) {
      REQUIRE(std::abs(out[i]) == Approx(0));
    }
  }
}

} // namespace mps