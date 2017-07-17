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
#include "binner.h"
#include "base.h"
#include "window.h"

namespace mps {

namespace {

Int GetTau(Int q, Int bins, Int idx) {
  if (idx == 0) {
    return q;
  }
  if (idx & 1) {
    return q + bins * (1 << ((idx - 1) / 2));
  }
  return q - bins * (1 << (idx / 2 - 1));
}

} // namespace

BinInTime::BinInTime(const Window &win, Int bits) : win_(win), bits_(bits) {
  plan_.reset(new FFTPlan(win.bins(), FFTW_FORWARD));
}

BinInTime *BinInTime::Create(int binner_type, const Window &win, Int bits) {
  if (binner_type == 0) {
    return new BinInTimeV0(win, bits);
  }
  if (binner_type == 1) {
    return new BinInTimeV1(win, bits);
  }
  if (binner_type == 2) {
    return new BinInTimeV2(win, bits);
  }
  LOG(FATAL) << "Unknown binner_type=" << binner_type;
  return nullptr;
}

BinInTimeV0::BinInTimeV0(const Window &win, Int bits)
    : BinInTime(win, bits), scratch_(win.bins()) {}

void BinInTimeV0::Run(const CplexArray &x, const Transform &tf, Int q,
                      CplexMatrix *out) {
  CHECK_EQ(out->rows(), 1 + 2 * bits_);
  const Int bins = win_.bins();
  const Int n = win_.n();

  const double delta = -0.5 / double(bins);
  DCHECK_EQ(n, x.size());

  const Int p = win_.p();
  const Int p2 = (p - 1) / 2;

  for (Int u = 0; u < out->rows(); ++u) {
    const Int tau = GetTau(q, bins, u);
    scratch_.Clear();
    for (Int i = 0; i < p; ++i) {
      const Int t = i <= p2 ? i : i - p;
      const double wt = win_.wt(i <= p2 ? i : p - i);
      const Int j = PosMod(tf.a * (t + tau) + tf.c, n);
      const Int k = Mod(tf.b * (t + tau), n);
      // Do mods for better accuracy. Note fmod can be negative, but it is ok.
      const double freq =
          (double(k) / double(n) + std::fmod(delta * double(t), 1.0));
      scratch_[i % bins] += (x[j] * Sinusoid(freq)) * wt;
    }
    // Do B-point FFT.
    CplexArray &v = (*out)[u];
    DCHECK_EQ(bins, v.size());
    plan_->Run(scratch_, &v);
  }
}

BinInTimeV1::BinInTimeV1(const Window &win, Int bits)
    : BinInTime(win, bits), scratch_(win.bins()), scratch2_(win.bins()) {}

inline void BinInTimeV1Helper(Int i, Int t, double wt, Int offset, Int q, Int n,
                              Int bins, Cplex bq, double delta,
                              const Transform &tf, const CplexArray &x,
                              CplexArray *scratch, CplexArray *scratch2) {
  const Int ts = t + offset;
  const Int j1 = PosMod(tf.a * (q + ts) + tf.c, n);
  const Int j2 = PosMod(tf.a * (q - ts) + tf.c, n);
  const Cplex x1 = x[j1];
  const Cplex x2 = x[j2];

  const Int k = Mod(tf.b * ts, n);
  // Do mods for better accuracy. Note fmod can be negative, but it is ok.
  const double freq = double(k) / double(n) + std::fmod(delta * double(t), 1.0);
  // Only one Sinusoid per iteration.
  const Cplex factor = Sinusoid(freq);

  const Cplex z1 = x1 * bq * factor * wt;
  const Cplex z2 = std::conj(x2) * std::conj(bq) * factor * wt;

  (*scratch)[i % bins] += 0.5 * (z1 + z2);
  (*scratch2)[i % bins] += RotateBackward(0.5 * (z1 - z2));
}

void BinInTimeV1::Run(const CplexArray &x, const Transform &tf, Int q,
                      CplexMatrix *out) {
  CHECK_EQ(out->rows(), 1 + 2 * bits_);
  const Int bins = win_.bins();
  const Int n = win_.n();

  const double delta = -0.5 / double(bins);
  DCHECK_EQ(n, x.size());

  const Int p = win_.p();
  const Int p2 = (p - 1) / 2;

  // Take care of tau=q first.
  scratch_.Clear();
  for (Int i = 0; i < p; ++i) {
    const Int t = i <= p2 ? i : i - p;
    const double wt = win_.wt(i <= p2 ? i : p - i);
    const Int j = PosMod(tf.a * (t + q) + tf.c, n);
    const Int k = Mod(tf.b * (t + q), n);
    const double freq =
        double(k) / double(n) + std::fmod(delta * double(t), 1.0);
    scratch_[i % bins] += (x[j] * Sinusoid(freq)) * wt;
  }
  plan_->Run(scratch_, &(*out)[0]);

  // Take care of offsets from q.
  const Cplex bq = Sinusoid(double(Mod(tf.b * q, n)) / double(n));

  for (Int bit = 0; bit < bits_; ++bit) {
    const Int offset = bins * (1 << bit);
    scratch_.Clear();
    scratch2_.Clear();
    for (Int i = 0; i <= p2; ++i) {
      const Int t = i;
      const double wt = win_.wt(i);
      BinInTimeV1Helper(i, t, wt, offset, q, n, bins, bq, delta, tf, x,
                        &scratch_, &scratch2_);
    }
    for (Int i = p2 + 1; i < p; ++i) {
      const Int t = i - p;
      const double wt = win_.wt(p - i);
      BinInTimeV1Helper(i, t, wt, offset, q, n, bins, bq, delta, tf, x,
                        &scratch_, &scratch2_);
    }

    // Do B-point FFT.
    plan_->Run(scratch_, &(*out)[2 * bit + 1]);
    plan_->Run(scratch2_, &(*out)[2 * bit + 2]);

    // Combine the bin coefficients.
    for (Int bin = 0; bin < bins; ++bin) {
      const Cplex coef_r = (*out)[2 * bit + 1][bin];
      const Cplex coef_i = (*out)[2 * bit + 2][bin];
      (*out)[2 * bit + 1][bin] = coef_r + RotateForward(coef_i);
      (*out)[2 * bit + 2][bin] =
          std::conj(coef_r) + RotateForward(std::conj(coef_i));
    }
  }
}

BinInTimeV2::BinInTimeV2(const Window &win, Int bits)
    : BinInTime(win, bits), scratch_(win.p()), scratch2_(win.p()),
      idx_(win.p()), idx2_(win.p()) {}

inline void BinInTimeV2Helper1(Int i, Int t, const Transform &tf, Int offset,
                               Int q, Int n, Int *idx, Int *idx2) {
  const Int ts = t + offset;
  idx[i] = PosMod(tf.a * (q + ts) + tf.c, n);
  idx2[i] = PosMod(tf.a * (q - ts) + tf.c, n);
}

inline void BinInTimeV2Helper2(Int i, Cplex bq, const CplexArray &x, Int *idx,
                               Int *idx2, Cplex *out, Cplex *out2) {
  const Cplex x1 = x[idx[i]] * bq;
  const Cplex x2 = std::conj(x[idx2[i]] * bq);
  out[i] = 0.5 * (x1 + x2);
  out2[i] = RotateBackward(0.5 * (x1 - x2));
}

inline void BinInTimeV2Helper3(Int i, Int t, const Transform &tf, Int offset,
                               double wt, Int n, double delta, Cplex *out,
                               Cplex *out2) {
  const Int ts = t + offset;
  const Int k = Mod(tf.b * ts, n);
  // Do mods for better accuracy. Note fmod can be negative, but it is ok.
  const double freq = double(k) / double(n) + std::fmod(delta * double(t), 1.0);
  const Cplex factor = Sinusoid(freq) * wt;
  out[i] *= factor;
  out2[i] *= factor;
}

// #pragma GCC push_options
// #pragma GCC optimize ("unroll-loops")
// #pragma GCC pop_options

void BinInTimeV2::Run(const CplexArray &x, const Transform &tf, Int q,
                      CplexMatrix *out) {
  CHECK_EQ(out->rows(), 1 + 2 * bits_);
  const Int bins = win_.bins();
  const Int n = win_.n();

  const double delta = -0.5 / double(bins);
  DCHECK_EQ(n, x.size());

  const Int p = win_.p();
  const Int p2 = (p - 1) / 2;

  // Take care of tau=q first. Treat this case simply.
  scratch_.Clear();
  for (Int i = 0; i < p; ++i) {
    const Int t = i <= p2 ? i : i - p;
    const double wt = win_.wt(i <= p2 ? i : p - i);
    const Int j = PosMod(tf.a * (t + q) + tf.c, n);
    const Int k = Mod(tf.b * (t + q), n);
    const double freq =
        double(k) / double(n) + std::fmod(delta * double(t), 1.0);
    scratch_[i % bins] += (x[j] * Sinusoid(freq)) * wt;
  }
  plan_->Run(scratch_, &(*out)[0]);

  // Take care of offsets from q.
  const Cplex bq = Sinusoid(double(Mod(tf.b * q, n)) / double(n));

  for (Int bit = 0; bit < bits_; ++bit) {
    const Int offset = bins * (1 << bit);

    for (Int i = 0; i <= p2; ++i) {
      BinInTimeV2Helper1(i, i, tf, offset, q, n, idx_.data(), idx2_.data());
    }
    for (Int i = p2 + 1; i < p; ++i) {
      BinInTimeV2Helper1(i, i - p, tf, offset, q, n, idx_.data(), idx2_.data());
    }

    for (Int i = 0; i < p; ++i) {
      BinInTimeV2Helper2(i, bq, x, idx_.data(), idx2_.data(), scratch_.data(),
                         scratch2_.data());
    }

    for (Int i = 0; i <= p2; ++i) {
      BinInTimeV2Helper3(i, i, tf, offset, win_.wt(i), n, delta,
                         scratch_.data(), scratch2_.data());
    }
    for (Int i = p2 + 1; i < p; ++i) {
      BinInTimeV2Helper3(i, i - p, tf, offset, win_.wt(p - i), n, delta,
                         scratch_.data(), scratch2_.data());
    }

    const Int folds = p / bins;
    for (Int i = 1; i < folds; ++i) {
      for (Int j = 0; j < bins; ++j) {
        const Int k = i * bins + j;
        scratch_[j] += scratch_[k];
        scratch2_[j] += scratch2_[k];
      }
    }

    // Do B-point FFT.
    plan_->Run(scratch_, &(*out)[2 * bit + 1]);
    plan_->Run(scratch2_, &(*out)[2 * bit + 2]);

    // Combine the bin coefficients.
    for (Int bin = 0; bin < bins; ++bin) {
      const Cplex coef_r = (*out)[2 * bit + 1][bin];
      const Cplex coef_i = (*out)[2 * bit + 2][bin];
      (*out)[2 * bit + 1][bin] = coef_r + RotateForward(coef_i);
      (*out)[2 * bit + 2][bin] =
          std::conj(coef_r) + RotateForward(std::conj(coef_i));
    }
  }
}

BinInFreq::BinInFreq(const Window &win, Int bits) : win_(win), bits_(bits) {}

BinInFreq *BinInFreq::Create(int binner_type, const Window &win, Int bits) {
  if (binner_type == 0) {
    return new BinInFreqV0(win, bits);
  }
  if (binner_type == 1) {
    return new BinInFreqV1(win, bits);
  }
  LOG(FATAL) << "Unknown binner_type=" << binner_type;
  return nullptr;
}

BinInFreqV0::BinInFreqV0(const Window &win, Int bits) : BinInFreq(win, bits) {}

void BinInFreqV0::Run(const ModeMap &mm, const Transform &tf, Int q,
                      CplexMatrix *out) {
  CHECK_EQ(out->rows(), 1 + 2 * bits_);
  const Int bins = win_.bins();
  const Int n = win_.n();

  for (const auto &kv : mm) {
    const Int k = kv.first;
    const Int l = PosMod(tf.a * k + tf.b, n); // 0 to n-1.
    const Int bin = Int(l * bins / n);
    const double xi =
        (double(bin) + 0.5) / double(bins) - double(l) / double(n);
    const double wf = win_.SampleInFreq(xi);
    const Cplex base = kv.second * wf;
    for (Int u = 0; u < out->rows(); ++u) {
      const Int tau = GetTau(q, bins, u);
      const Int s = Mod(tf.c * k + l * tau, n);
      (*out)[u][bin] -= base * Sinusoid(double(s) / double(n));
    }
  }
}

BinInFreqV1::BinInFreqV1(const Window &win, Int bits) : BinInFreq(win, bits) {}

void BinInFreqV1::Run(const ModeMap &mm, const Transform &tf, Int q,
                      CplexMatrix *out) {
  CHECK_EQ(out->rows(), 1 + 2 * bits_);
  const Int bins = win_.bins();
  const Int n = win_.n();

  for (const auto &kv : mm) {
    const Int k0 = kv.first;
    const Int k1 = PosMod(tf.a * k0 + tf.b, n); // 0 to n-1.
    const Int bin = Int(k1 * bins / n);
    const double xi =
        (double(bin) + 0.5) / double(bins) - double(k1) / double(n);
    const double wf = win_.SampleInFreq(xi);

    const double freq = double(Mod(q * k1 + tf.c * k0, n)) / double(n);
    const Cplex base = wf * kv.second * Sinusoid(freq);
    (*out)[0][bin] -= base;

    for (Int bit = 0; bit < bits_; ++bit) {
      const Int offset = bins * (1 << bit);
      const Cplex factor = Sinusoid(double(Mod(offset * k1, n)) / double(n));
      (*out)[bit * 2 + 1][bin] -= base * factor;
      (*out)[bit * 2 + 2][bin] -= base * std::conj(factor); // Divide by factor.
    }
  }
}

} // namespace mps