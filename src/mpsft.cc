#include "mpsft.h"
#include "freqid.h"

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
           Real threshold, ModeMap* mm) {
  CHECK_EQ(1, trials % 2) << "Odd number of trials expected";

  const Int n = x.size();
  Window win(n, bins, delta);
  const Int bits = NumBits(n, bins);

  Transform tf(n);
  FFTPlan plan(bins, FFTW_FORWARD);
  CplexArray scratch(bins);

  TauSet taus;
  for (Int bit = 0; bit < bits; ++bit) {
    taus.list_s.push_back((1 << bit) * bins);
  }
  const Int m = taus.size();

  vector<unique_ptr<CplexMatrix>> bin_coefs(trials);
  // Repeat "trials" number of times.
  for (Int trial = 0; trial < trials; ++trial) {
    taus.q = RandomInt() % n;
    bin_coefs[trial].reset(new CplexMatrix(m, bins));
    // BinInTime will produce "bins" number of coefficients for each tau in
    // taus.
    BinInTime(win, tf, taus, x, &plan, bin_coefs[trial].get(), &scratch);
  }

  Real sigma[2];
  for (Int b = 0; b < bins; ++b) {
    Cplex sum = 0;
    for (Int trial = 0; trial < trials; ++trial) {
      sum += (*bin_coefs[trial])[0][b];
    }
    sum /= trials;
    const Real bin_energy = AbsSq(sum);
    if (bin_energy < threshold * threshold) {
      continue;
    }

    Int xi_sum = 0;
    for (Int bit = 0; bit < bits; ++bit) {
      Int count = 0;
      for (Int trial = 0; trial < trials; ++trial) {
        CplexMatrix &a = *bin_coefs[trial]; // Bin coefficients.
        const Cplex u1 = a[2 * bit + 1][b];
        const Cplex u2 = a[2 * bit + 2][b];
        const Cplex u0 = a[0][b];
        if (MatPencil(u0, u1, u2, sigma)) {
          ++count;
        }
      }
      LOG(INFO) << "bit=" << bit << " count=" << count;
      xi_sum <<= 1;
      if (count > trials / 2) {
        ++xi_sum;
      }
    }
    // xi is between 0 and 1.
    // [0, 1] is divided into 2^bits minibins.
    // xi is the center of one of these minibins.
    const Real xi = Real(2 * xi_sum + 1) / Real(1 << (bits + 1));
    const Int k1 = std::round(Real(n) * (Real(b) + xi) / Real(bins));
    const Int k0 = PosMod(tf.a_inv * (k1 - tf.b), n);

    // Estimate the coefficient.
  }
}

} // namespace mps