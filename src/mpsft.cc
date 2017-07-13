#include "mpsft.h"

namespace mps {

namespace {

// Returns 1 + floor(log2(n / bins))
Int NumBits(Int n, Int bins) {
  if (n < bins) {
    return 1;
  }
  Int out = 0;
  // n=2*bins should return 2.
  // n=2*bins+1 should return 2.
  while (n >= bins) {
    n >>= 1;
    ++out;
  }
  return out;
}

} // namespace

void mpsft(const CplexArray &x, Int bins, Real delta, Int trials,
           Real threshold) {
  CHECK_EQ(1, trials % 2) << "Odd number of trials expected";

  const Int n = x.size();
  Window win(n, bins, delta);
  const Int bits = NumBits(n, bins);

  Transform tf(n);
  FFTPlan plan(bins, FFTW_FORWARD);
  CplexArray scratch(bins);

  TauSet taus;
  for (Int i = 0; i < bits; ++i) {
    taus.list_s.push_back((1 << i) * bins);
  }
  const Int m = taus.size();

  vector<unique_ptr<CplexMatrix>> bin_coefs(trials);
  for (Int i = 0; i < trials; ++i) {
    taus.q = RandomInt() % n;
    bin_coefs[i].reset(new CplexMatrix(m, bins));
    BinInTime(win, tf, taus, x, &plan, bin_coefs[i].get(), &scratch);
  }

  for (Int b = 0; b < bins; ++b) {
    Cplex sum = 0;
    for (Int j = 0; j < trials; ++j) {
      sum += (*bin_coefs[j])[0][b];
    }
    sum /= trials;
    const Real bin_energy = AbsSq(sum);
    if (bin_energy < threshold * threshold) {
      continue;
    }
    LOG(INFO) << "b=" << b;
  }
}

} // namespace mps