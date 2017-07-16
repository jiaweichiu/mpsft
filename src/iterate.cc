/*
 * Copyright (c) 2017 Jiawei Chiu
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301 USA
 *
 */
#include "iterate.h"
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

std::pair<double, bool>
IdentifyFreq(const IterateOptions &opt, Int b, Int bits,
             const vector<std::unique_ptr<CplexMatrix>> &bin_coefs) {
  const Int trials = opt.trials;
  double sigma[2];

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
      // 0.5 factor is due to matrix being 2 by 2.
      if (opt.sv_threshold > 0 && sigma[1] * 0.5 > opt.sv_threshold) {
        return std::make_pair(0, false);
      }
    }
    // LOG(INFO) << "bit=" << bit << " count=" << count;
    xi_sum <<= 1;
    if (count > trials / 2) {
      ++xi_sum;
    }
  }
  // xi is between 0 and 1.
  // [0, 1] is divided into 2^bits minibins.
  // xi is the center of one of these minibins.
  return std::make_pair(double(2 * xi_sum + 1) / double(1 << (bits + 1)), true);
}

} // namespace

bool Iterate(const CplexArray &x, const IterateOptions &opt, ModeMap *mm) {
  bool found = false; // Whether we found anything new.
  const Int trials = opt.trials;
  const Int bins = opt.bins;
  CHECK_GT(opt.window_threshold, 0);

  CHECK_EQ(1, trials % 2) << "Odd number of trials expected";

  const Int n = x.size();
  Window win(n, bins, opt.window_delta);
  const Int bits = NumBits(n, bins);

  Transform tf(n);
  BinnerSimple binner(win, bits);

  vector<double> list_q;
  vector<std::unique_ptr<CplexMatrix>> bin_coefs(trials);

  // Repeat "trials" number of times.
  for (Int trial = 0; trial < trials; ++trial) {
    const Int q = RandomInt() % n;
    list_q.push_back(q);

    bin_coefs[trial].reset(new CplexMatrix(1 + 2 * bits, bins));
    // BinInTime will produce "bins" number of coefficients for each tau.
    CplexMatrix *a = bin_coefs[trial].get();
    binner.BinInTime(x, tf, q, a);
    binner.BinInFreq(*mm, tf, q, a);
  }

  for (Int b = 0; b < bins; ++b) {
    if (opt.bin_threshold > 0) {
      double bin_energy = 0;
      for (Int trial = 0; trial < trials; ++trial) {
        bin_energy += AbsSq((*bin_coefs[trial])[0][b]);
      }
      bin_energy /= trials;
      if (bin_energy < opt.bin_threshold * opt.bin_threshold) {
        continue;
      }
    }

    // res.first is xi, in [0, 1].
    auto res = IdentifyFreq(opt, b, bits, bin_coefs);
    if (!res.second) {
      continue;
    }
    const double xi1 = res.first;

    // k1 is mode location after random permutation.
    const Int k1 = std::round(double(n) * (double(b) + xi1) / double(bins));
    DCHECK_GE(k1, 0);
    DCHECK_LT(k1, n);

    // Check window in frequency / attentuation factor. If too small, reject.
    const double xi0 =
        double(k1) / double(n) - (double(b) + 0.5) / double(bins);
    const double wf = win.SampleInFreq(xi0);
    if (std::abs(wf) < opt.window_threshold) {
      continue;
    }

    // Estimate coefficient in transformed signal.
    Cplex coef_sum = 0;
    for (Int trial = 0; trial < trials; ++trial) {
      const CplexMatrix &a = *bin_coefs[trial];
      const double freq = -double(Mod(list_q[trial] * k1, n)) / double(n);
      const Cplex factor = Sinusoid(freq);
      coef_sum += a[0][b] * factor;

      for (Int bit = 0; bit < bits; ++bit) {
        // Try to do only one Sinusoid here instead of two, using symmetry.
        const double freq2 =
            -double(Mod(bins * (1 << bit) * k1, n)) / double(n);
        const Cplex factor2 = Sinusoid(freq2);
        const Cplex f1 = factor * factor2;
        const Cplex f2 = factor * std::conj(factor2); // Divide by factor2.
        coef_sum += a[2 * bit + 1][b] * f1;
        coef_sum += a[2 * bit + 2][b] * f2;
      }
    }
    Cplex coef = coef_sum / double(trials * (1 + 2 * bits));

    // Undo the transform.
    // k0 is original mode location.
    const Int k0 = PosMod(tf.a_inv * (k1 - tf.b), n);
    coef *= Sinusoid(-double(Mod(tf.c * k0, n)) / double(n));
    coef /= wf;

    (*mm)[k0] += coef;
    found = true;
  }
  return found;
}

} // namespace mps