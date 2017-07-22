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
#include "base.h"
#include "freqid.h"
#include "sincos.h"

namespace mps {

namespace {

// Returns 1 + floor(log2(n / bins))
int32_t NumBits(int32_t n, int32_t bins) {
  if (n < bins) {
    return 1;
  }
  int32_t out = 0;
  // n=2*bins should return 2.
  // n=2*bins+1 should return 2.
  while (n >= bins) {
    n >>= 1;
    ++out;
  }
  return out;
}

std::pair<double, bool>
IdentifyFreq(const IterateOptions &opt, int32_t b, int32_t bits,
             const vector<std::unique_ptr<CplexMatrix>> &bin_coefs) {
  const int32_t trials = opt.trials;
  double sigma[2];

  int32_t xi_sum = 0;
  for (int32_t bit = 0; bit < bits; ++bit) {
    int32_t count = 0;
    for (int32_t trial = 0; trial < trials; ++trial) {
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
  const int32_t trials = opt.trials;
  const int32_t bins = opt.bins;
  CHECK_GT(opt.window_threshold, 0);

  CHECK_EQ(1, trials % 2) << "Odd number of trials expected";

  const int32_t n = x.size();
  Window win(n, bins, opt.window_delta);
  const int32_t bits = NumBits(n, bins);

  Transform tf(n);
  std::unique_ptr<BinInTime> bin_in_time(
      BinInTime::Create(opt.bin_in_time_type, win, bits));
  std::unique_ptr<BinInFreq> bin_in_freq(
      BinInFreq::Create(opt.bin_in_freq_type, win, bits));

  vector<double> list_q;
  vector<std::unique_ptr<CplexMatrix>> bin_coefs(trials);

  // Repeat "trials" number of times.
  for (int32_t trial = 0; trial < trials; ++trial) {
    const int32_t q = RandomInt32() % n;
    list_q.push_back(q);

    bin_coefs[trial].reset(new CplexMatrix(1 + 2 * bits, bins));
    // BinInTime will produce "bins" number of coefficients for each tau.
    CplexMatrix *a = bin_coefs[trial].get();
    bin_in_time->Run(x, tf, q, a);
    bin_in_freq->Run(*mm, tf, q, a);
  }

  for (int32_t b = 0; b < bins; ++b) {
    if (opt.bin_threshold > 0) {
      double bin_energy = 0;
      for (int32_t trial = 0; trial < trials; ++trial) {
        const Cplex coef = (*bin_coefs[trial])[0][b];
        bin_energy += AbsSq(RE(coef), IM(coef));
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
    const int32_t k1 = std::round(double(n) * (double(b) + xi1) / double(bins));
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
    for (int32_t trial = 0; trial < trials; ++trial) {
      const CplexMatrix &a = *bin_coefs[trial];
      const double freq = -double(MulMod(list_q[trial], k1, n)) / double(n);
      const Cplex factor = Sinusoid(freq);
      coef_sum += a[0][b] * factor;

      for (int32_t bit = 0; bit < bits; ++bit) {
        // Try to do only one Sinusoid here instead of two, using symmetry.
        const double freq2 =
            -double(MulMod(bins * (1 << bit), k1, n)) / double(n);
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
    const int32_t k0 = MulPosMod(tf.a_inv, int64_t(k1) - int64_t(tf.b), n);
    coef *= Sinusoid(-double(MulMod(tf.c, k0, n)) / double(n));
    coef /= wf;

    (*mm)[k0] += coef;
    found = true;
  }
  return found;
}

} // namespace mps