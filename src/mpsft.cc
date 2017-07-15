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

void Iterate(const CplexArray &x, const IterateOptions &opt, ModeMap *mm) {
  const Int trials = opt.trials;
  const Int bins = opt.bins;

  CHECK_EQ(1, trials % 2) << "Odd number of trials expected";

  const Int n = x.size();
  Window win(n, bins, opt.window_delta);
  const Int bits = NumBits(n, bins);

  Transform tf(n);
  FFTPlan plan(bins, FFTW_FORWARD);
  CplexArray scratch(bins);

  vector<Real> list_q;
  TauSet taus;
  for (Int bit = 0; bit < bits; ++bit) {
    taus.list_s.push_back((1 << bit) * bins);
  }

  vector<unique_ptr<CplexMatrix>> bin_coefs(trials);
  // Repeat "trials" number of times.
  for (Int trial = 0; trial < trials; ++trial) {
    const Int q = RandomInt() % n;
    list_q.push_back(q);
    taus.q = q;
    bin_coefs[trial].reset(new CplexMatrix(taus.size(), bins));
    // BinInTime will produce "bins" number of coefficients for each tau in
    // taus.
    CplexMatrix *a = bin_coefs[trial].get();
    BinInTime(win, tf, taus, x, &plan, a, &scratch);
    BinInFreq(win, tf, taus, *mm, a);
  }

  Real sigma[2];
  for (Int b = 0; b < bins; ++b) {
    Real bin_energy = 0;
    for (Int trial = 0; trial < trials; ++trial) {
      bin_energy += AbsSq((*bin_coefs[trial])[0][b]);
    }
    bin_energy /= trials;
    if (bin_energy < opt.bin_threshold * opt.bin_threshold) {
      continue;
    }

    // LOG(INFO) << "b=" << b << " energy=" << std::sqrt(bin_energy);

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
      // LOG(INFO) << "bit=" << bit << " count=" << count;
      xi_sum <<= 1;
      if (count > trials / 2) {
        ++xi_sum;
      }
    }
    // xi1 is between 0 and 1.
    // [0, 1] is divided into 2^bits minibins.
    // xi1 is the center of one of these minibins.
    const Real xi1 = Real(2 * xi_sum + 1) / Real(1 << (bits + 1));

    // k1 is mode location after random permutation.
    const Int k1 = std::round(Real(n) * (Real(b) + xi1) / Real(bins));
    DCHECK_GE(k1, 0);
    DCHECK_LT(k1, n);

    // Check window in frequency / attentuation factor. If too small, reject.
    const Real xi = Real(k1) / Real(n) - (Real(b) + 0.5) / Real(bins);
    const Real wf = win.SampleInFreq(xi);

    if (std::abs(wf) < opt.window_threshold) {
      continue;
    }

    // Estimate coefficient in transformed signal.
    Cplex coef_sum = 0;
    for (Int trial = 0; trial < trials; ++trial) {
      const CplexMatrix &a = *bin_coefs[trial];

      const Real angle =
          -2.0 * M_PI * Real(Mod(Long(list_q[trial]) * Long(k1), n)) / Real(n);
      const Cplex factor = Sinusoid(angle);
      coef_sum += a[0][b] * factor;

      for (Int bit = 0; bit < bits; ++bit) {
        // Try to do only one Sinusoid here instead of two, using symmetry.
        const Real angle2 = -2.0 * M_PI *
                            Real(Mod(Long(taus.list_s[bit]) * Long(k1), n)) /
                            Real(n);
        const Cplex factor2 = Sinusoid(angle2);
        const Cplex f1 = factor * factor2;
        const Cplex f2 = factor * std::conj(factor2); // Divide by factor2.
        coef_sum += a[2 * bit + 1][b] * f1;
        coef_sum += a[2 * bit + 2][b] * f2;
      }
    }
    Cplex coef = coef_sum / Real(trials * taus.size());

    // Undo the transform.
    // k0 is original mode location.
    const Int k0 = PosMod(Long(tf.a_inv) * (Long(k1) - Long(tf.b)), n);
    coef *=
        Sinusoid(-2.0 * M_PI * Real(Mod(Long(tf.c) * Long(k0), n)) / Real(n));
    coef /= wf;

    (*mm)[k0] += coef;
  }
}

} // namespace mps