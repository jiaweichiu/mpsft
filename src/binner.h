#pragma once

#include "base.h"
#include "window.h"

namespace mps {

struct TauSet {
  // tau's are: q, q +/- s for s in list_s.
  Int q;
  vector<Int> list_s;

  inline Int size() const { return 1 + 2 * list_s.size(); }

  // Tau(0) = q; Tau(1) = q+s[0]; Tau(2) = q-s[0]; ...
  inline Int value(Int idx) const {
    if (idx == 0) {
      return q;
    }
    if (idx & 1) { // Odd.
      return q + list_s[(idx - 1) / 2];
    }
    return q - list_s[idx / 2 - 1];
  }
};

/*
1) Assume plan is forward FFT and has size equal to number of bins.
2) We will clear out1, use it as scratch and output in out2.
3) Number of tau's is equal to 1+2*len(list_s).
4) len(out) is equal to number of tau's. Each element is a CplexArray of size B.
5) scratch is of size B where B is number of bins.
6) out[0] -> q, out[1] -> q+s[0], out[2] -> q-s[0], out[3] -> q+s[1], ...
*/
void BinInTime(const Window &win, const Transform &tf, const TauSet &taus,
               const CplexArray &x, FFTPlan *plan, CplexMatrix *out,
               CplexArray *scratch);

// Subtract results from out. Similar convention as BinInTime.
void BinInFreq(const Window &win, const Transform &tf, const TauSet &taus,
               const CplexArray &coef, const vector<Int> &loc, CplexMatrix *out);

} // namespace mps