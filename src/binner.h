#pragma once

#include "base.h"
#include "window.h"

namespace mps {

/*struct TauSet {
  // tau's are: q, q +/- s for s in list_s.
  Int q;
  vector<Int> list_s;

  inline Int size() const { return 1 + 2 * list_s.size(); }

  // Tau(0) = q; Tau(1) = q+s[0]; Tau(2) = q-s[0]; ...
  inline Int value(Int idx) const;
};*/

/*class TauSet {
public:
  TauSet(Int q, Int bins, Int bits);
  inline Int size() const { return 1 + 2 * bits_; }
  inline Int tau(Int idx) const;
  inline Int q() const { return q_; }
  inline Int offset(Int i) const { return bins_ * (1 << i); }

private:
  Int q_;
  Int bins_;
  Int bits_;
};*/

/*
1) Assume plan is forward FFT and has size equal to number of bins.
2) We will clear out1, use it as scratch and output in out2.
3) Number of tau's is equal to 1+2*len(list_s).
4) len(out) is equal to number of tau's. Each element is a CplexArray of size B.
5) scratch is of size B where B is number of bins.
6) out[0] -> q, out[1] -> q+s[0], out[2] -> q-s[0], out[3] -> q+s[1], ...
*/
class Binner {
public:
  Binner(const Window &win, const Transform &tf, Int bits);

  void BinInTime(const CplexArray &x, Int q, CplexMatrix *out);
  void BinInFreq(const ModeMap &mm, Int q, CplexMatrix *out);

private:
  const Window &win_;
  const Transform &tf_;
  Int bits_;
  unique_ptr<FFTPlan> plan_;
  unique_ptr<CplexArray> scratch_;
};

} // namespace mps