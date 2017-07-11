#pragma once

#include "base.h"
#include "window.h"

namespace mps {

struct TauSet {
  // tau's are: q, q +/- s for s in list_s.
  Int q;
  vector<Int> list_s;
};

// Assume plan is forward FFT and has size equal to number of bins.
// We will clear out1, use it as scratch and output in out2.
void BinInTime(const Window &win, const Transform &tf, const TauSet &taus,
               const CplexArray &x, FFTPlan *plan, CplexArray *out1,
               CplexArray *out2);

// Will subtract from out.
void BinInFreq(const Window &win, const Transform &tf, const TauSet &taus,
               const CplexArray &coef, const vector<Int> &loc, CplexArray *out);

} // namespace mps