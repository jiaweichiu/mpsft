#pragma once

#include "base.h"
#include "window.h"

namespace mps {

class Binner {
public:
  Binner(const Window &win, const Transform &tf, Int bits);

  // tau = q +/- (1 << b) where 0 <= b < bits.
  // Num of rows of "out" is 2*bits+1.
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