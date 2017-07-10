#pragma once

#include "base.h"
#include "perm.h"

namespace mps {

struct BinnerOpt {
  Int n;    // Assumed to be prime.
  Int bins; // Number of bins.
  Real delta;
};

struct TauSet {
  // tau's are: q, q +/- s for s in list_s.
  Int q;
  vector<Int> list_s;
};

void BinInTime(const Window &win, const Transform &tf, const TauSet &taus,
               const vector<Cplex> &x);

} // namespace mps