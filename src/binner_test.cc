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
  for (Int i = 0; i < coef.Size(); ++i) {
    for (Int t = 0; t < n; ++t) {
      x[t] += coef[i] * Sinusoid((2.0 * M_PI) * Real(t) / Real(n));
    }
  }

  // Prepare binning.
  Transform tf(n, 3, 847, 45);
  TauSet taus;
  taus.q = 106;

  // BinInTime.
  CplexArray out1(bins);
  CplexArray out2(bins);
  FFTPlan plan(bins, -1);
  
  BinInTime(win, tf, taus, x, &plan, &out1, &out2);
  for (Int i = 0; i < bins; ++i) {
    LOG(INFO) << out2[i];
  }

  // BinInFreq.
  CplexArray out(bins);
  out.Clear();
  BinInFreq(win, tf, taus, coef, loc, &out); // Subtract.
  for (Int i = 0; i < bins; ++i) {
    LOG(INFO) << out[i];
  }

  REQUIRE(RE(out[0]) == Approx(-1.12652));
  REQUIRE(IM(out[0]) == Approx(0.108838));
  for (Int i = 1; i < bins; ++i) {
    REQUIRE(std::abs(out[i]) == Approx(0));
  }
}

} // namespace mps