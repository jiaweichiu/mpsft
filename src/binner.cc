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
#include <boost/align/aligned_allocator.hpp>
#include <boost/simd/arithmetic.hpp>
#include <boost/simd/function/aligned_store.hpp>
#include <boost/simd/function/load.hpp>
#include <boost/simd/function/sin.hpp>
#include <boost/simd/function/sincos.hpp>
#include <boost/simd/function/store.hpp>
#include <boost/simd/memory/allocator.hpp>
#include <boost/simd/pack.hpp>

#include "base.h"
#include "binner.h"
#include "integer.h"
#include "sincos.h"
#include "window.h"

namespace ba = boost::alignment;
namespace bs = boost::simd;

namespace mps {

namespace {

int32_t GetTau(int32_t q, int32_t bins, int32_t idx) {
  if (idx == 0) {
    return q;
  }
  if (idx & 1) {
    return q + bins * (1 << ((idx - 1) / 2));
  }
  return q - bins * (1 << (idx / 2 - 1));
}

} // namespace

BinInTime::BinInTime(const Window &win, int32_t bits) : win_(win), bits_(bits) {
  plan_.reset(new FFTPlan(win.bins(), FFTW_FORWARD));
}

BinInTime *BinInTime::Create(int binner_type, const Window &win, int32_t bits) {
  if (binner_type == 0) {
    return new BinInTimeV0(win, bits);
  }
  if (binner_type == 1) {
    return new BinInTimeV1(win, bits);
  }
  if (binner_type == 2) {
    return new BinInTimeV2(win, bits);
  }
  if (binner_type == 3) {
    return new BinInTimeV3(win, bits);
  }
  LOG(FATAL) << "Unknown binner_type=" << binner_type;
  return nullptr;
}

BinInTimeV0::BinInTimeV0(const Window &win, int bits)
    : BinInTime(win, bits), scratch_(win.bins()) {}

void BinInTimeV0::Run(const CplexArray &x, const Transform &tf, int32_t q,
                      CplexMatrix *out) {
  CHECK_EQ(out->rows(), 1 + 2 * bits_);
  const int32_t bins = win_.bins();
  const int32_t n = win_.n();

  const double delta = -0.5 / double(bins);
  DCHECK_EQ(n, x.size());

  const int32_t p = win_.p();
  const int32_t p2 = (p - 1) / 2;

  for (int32_t u = 0; u < out->rows(); ++u) {
    const int32_t tau = GetTau(q, bins, u);
    scratch_.clear();
    for (int32_t i = 0; i < p; ++i) {
      const int32_t t = i <= p2 ? i : i - p;
      const double wt = win_.wt()[i];
      const int32_t j = MulAddPosMod(tf.a, int64_t(t) + int64_t(tau), tf.c, n);
      const int32_t k = MulMod(tf.b, int64_t(t) + int64_t(tau), n);
      // Do mods for better accuracy. Note fmod can be negative, but it is ok.
      const double freq =
          (double(k) / double(n) + PosModOne(delta * double(t)));
      scratch_[i % bins] += (x[j] * Sinusoid(freq)) * wt;
    }
    // Do B-point FFT.
    CplexArray &v = (*out)[u];
    DCHECK_EQ(bins, v.size());
    plan_->Run(scratch_, &v);
  }
}

BinInTimeV1::BinInTimeV1(const Window &win, int32_t bits)
    : BinInTime(win, bits), scratch_(win.bins()), scratch2_(win.bins()) {}

inline void BinInTimeV1Helper(int32_t i, int32_t t, double wt, int32_t offset,
                              int32_t q, int32_t n, int32_t bins, Cplex bq,
                              double delta, const Transform &tf,
                              const CplexArray &x, CplexArray *scratch,
                              CplexArray *scratch2) {
  const int32_t ts = t + offset;
  const int32_t j1 = MulAddPosMod(tf.a, int64_t(q) + int64_t(ts), tf.c, n);
  const int32_t j2 = MulAddPosMod(tf.a, int64_t(q) - int64_t(ts), tf.c, n);
  const Cplex x1 = x[j1];
  const Cplex x2 = x[j2];

  const int32_t k = MulMod(tf.b, ts, n);
  // Do mods for better accuracy. Note fmod can be negative, but it is ok.
  const double freq = double(k) / double(n) + std::fmod(delta * double(t), 1.0);
  // Only one Sinusoid per iteration.
  const Cplex factor = Sinusoid(freq);

  const Cplex z1 = x1 * bq * factor * wt;
  const Cplex z2 = std::conj(x2) * std::conj(bq) * factor * wt;

  (*scratch)[i % bins] += 0.5 * (z1 + z2);
  (*scratch2)[i % bins] += RotateBackward(0.5 * (z1 - z2));
}

void BinInTimeV1::Run(const CplexArray &x, const Transform &tf, int32_t q,
                      CplexMatrix *out) {
  CHECK_EQ(out->rows(), 1 + 2 * bits_);
  const int32_t bins = win_.bins();
  const int32_t n = win_.n();

  const double delta = -0.5 / double(bins);
  DCHECK_EQ(n, x.size());

  const int32_t p = win_.p();

  // Take care of tau=q first.
  scratch_.clear();
  for (int32_t i = 0; i < p; ++i) {
    const int32_t t = win_.t()[i];
    const double wt = win_.wt()[i];
    const int32_t j = MulAddPosMod(tf.a, int64_t(t) + int64_t(q), tf.c, n);
    const int32_t k = MulMod(tf.b, int64_t(t) + int64_t(q), n);
    const double freq =
        double(k) / double(n) + std::fmod(delta * double(t), 1.0);
    scratch_[i % bins] += (x[j] * Sinusoid(freq)) * wt;
  }
  plan_->Run(scratch_, &(*out)[0]);

  // Take care of offsets from q.
  const Cplex bq = Sinusoid(double(MulMod(tf.b, q, n)) / double(n));

  for (int32_t bit = 0; bit < bits_; ++bit) {
    const int32_t offset = bins * (1 << bit);
    scratch_.clear();
    scratch2_.clear();
    for (int32_t i = 0; i < p; ++i) {
      const int32_t t = win_.t()[i];
      const double wt = win_.wt()[i];
      BinInTimeV1Helper(i, t, wt, offset, q, n, bins, bq, delta, tf, x,
                        &scratch_, &scratch2_);
    }

    // Do B-point FFT.
    plan_->Run(scratch_, &(*out)[2 * bit + 1]);
    plan_->Run(scratch2_, &(*out)[2 * bit + 2]);

    // Combine the bin coefficients.
    for (int32_t bin = 0; bin < bins; ++bin) {
      const Cplex coef_r = (*out)[2 * bit + 1][bin];
      const Cplex coef_i = (*out)[2 * bit + 2][bin];
      (*out)[2 * bit + 1][bin] = coef_r + RotateForward(coef_i);
      (*out)[2 * bit + 2][bin] =
          std::conj(coef_r) + RotateForward(std::conj(coef_i));
    }
  }
}

BinInTimeV2::BinInTimeV2(const Window &win, int32_t bits)
    : BinInTime(win, bits), s1_(win.p()), s2_(win.p()), idx1_(win.p()),
      idx2_(win.p()), dbl_(win.p()) {}

void BinInTimeV2::Run(const CplexArray &x, const Transform &tf, int32_t q,
                      CplexMatrix *out) {
  int32_t *__restrict__ t = win_.t();
  double *__restrict__ wt = win_.wt();
  Cplex *__restrict__ s1 = s1_.data();
  Cplex *__restrict__ s2 = s2_.data();
  int32_t *__restrict__ idx1 = idx1_.data();
  int32_t *__restrict__ idx2 = idx2_.data();
  double *__restrict__ dbl = dbl_.data();

  CHECK_EQ(out->rows(), 1 + 2 * bits_);
  const int32_t bins = win_.bins();
  const int32_t n = win_.n();

  const double delta = -0.5 / double(bins);
  DCHECK_EQ(n, x.size());

  const int32_t p = win_.p();

  // Take care of tau=q first. TODO: Optimize this later.
  s1_.clear();
  for (int32_t i = 0; i < p; ++i) {
    const int32_t tt = t[i];
    const int32_t j = MulAddPosMod(tf.a, int64_t(tt) + int64_t(q), tf.c, n);
    const int32_t k = MulMod(tf.b, int64_t(tt) + int64_t(q), n);
    const double freq = double(k) / double(n) + PosModOne(delta * double(tt));
    s1[i % bins] += (x[j] * Sinusoid(freq)) * wt[i];
  }
  plan_->Run(s1_, &(*out)[0]);

  // Take care of offsets from q.
  const Cplex bq = Sinusoid(double(MulMod(tf.b, q, n)) / double(n));

  for (int32_t bit = 0; bit < bits_; ++bit) {
    const int32_t offset = bins * (1 << bit);

#pragma omp simd aligned(idx1, idx2, dbl, t : kAlign)
    // int64 operations are not vectorized in AVX2.
    for (int32_t i = 0; i < p; ++i) {
      const int64_t ts = int64_t(t[i]) + int64_t(offset);
      idx1[i] = MulAddPosMod(tf.a, int64_t(q) + int64_t(ts), tf.c, n);
      idx2[i] = MulAddPosMod(tf.a, int64_t(q) - int64_t(ts), tf.c, n);
      const int32_t k1 = MulMod(tf.b, ts, n);
      dbl[i] = double(k1) / double(n) + PosModOne(delta * double(t[i]));
    }

#pragma omp simd aligned(idx1, idx2, s1, s2 : kAlign)
    // This doesn't seem vectorized either. Splitting into components helps.
    for (int32_t i = 0; i < p; ++i) {
      const Cplex x1 = x[idx1[i]] * bq;
      const Cplex x2 = std::conj(x[idx2[i]] * bq);
      s1[i] = 0.5 * (x1 + x2);
      s2[i] = RotateBackward(0.5 * (x1 - x2));
    }

// CAUTION: Make sure the following is vectorized.
#pragma omp simd aligned(s1, s2, dbl, t, wt : kAlign)
    for (int32_t i = 0; i < p; ++i) {
      const double freq = PosModOne(dbl[i]);
      const double cc = wt[i] * CosTwoPiApprox(freq);
      const double ss = wt[i] * SinTwoPiApprox(freq);

      s1[i] *= Cplex(cc, ss);
      s2[i] *= Cplex(cc, ss);
    }

    const int32_t folds = p / bins;
    for (int32_t i = 1; i < folds; ++i) {
      for (int32_t j = 0; j < bins; ++j) {
        const int32_t k = i * bins + j;
        s1[j] += s1[k];
        s2[j] += s2[k];
      }
    }

    // Do B-point FFT.
    plan_->Run(s1_, &(*out)[2 * bit + 1]);
    plan_->Run(s2_, &(*out)[2 * bit + 2]);

    // Combine the bin coefficients.
    for (int32_t bin = 0; bin < bins; ++bin) {
      const Cplex coef_r = (*out)[2 * bit + 1][bin];
      const Cplex coef_i = (*out)[2 * bit + 2][bin];
      (*out)[2 * bit + 1][bin] = coef_r + RotateForward(coef_i);
      (*out)[2 * bit + 2][bin] =
          std::conj(coef_r) + RotateForward(std::conj(coef_i));
    }
  }
}

BinInTimeV3::BinInTimeV3(const Window &win, int32_t bits)
    : BinInTime(win, bits), s1_(win.p()), s2_(win.p()), idx1_(win.p()),
      idx2_(win.p()), d1_(win.p()), d2_(win.p()), d3_(win.p()) {}

void BinInTimeV3::Run(const CplexArray &x, const Transform &tf, int32_t q,
                      CplexMatrix *out) {
  int32_t *__restrict__ t = win_.t();
  double *__restrict__ wt = win_.wt();
  Cplex *__restrict__ s1 = s1_.data();
  Cplex *__restrict__ s2 = s2_.data();
  int32_t *__restrict__ idx1 = idx1_.data();
  int32_t *__restrict__ idx2 = idx2_.data();
  double *__restrict__ d1 = d1_.data();
  double *__restrict__ d2 = d2_.data();
  double *__restrict__ d3 = d3_.data();
  using pack_t = bs::pack<double>;
  const size_t pack_card = bs::cardinal_of<pack_t>();

  CHECK_EQ(out->rows(), 1 + 2 * bits_);
  const int32_t bins = win_.bins();
  const int32_t n = win_.n();

  const double delta = -0.5 / double(bins);
  DCHECK_EQ(n, x.size());

  const int32_t p = win_.p();

  // Take care of tau=q first. TODO: Optimize this later.
  s1_.clear();
  for (int32_t i = 0; i < p; ++i) {
    const int32_t tt = t[i];
    const int32_t j = MulAddPosMod(tf.a, int64_t(tt) + int64_t(q), tf.c, n);
    const int32_t k = MulMod(tf.b, int64_t(tt) + int64_t(q), n);
    const double freq = double(k) / double(n) + PosModOne(delta * double(tt));
    s1[i % bins] += (x[j] * Sinusoid(freq)) * wt[i];
  }
  plan_->Run(s1_, &(*out)[0]);

  // Take care of offsets from q.
  const Cplex bq = Sinusoid(double(MulMod(tf.b, q, n)) / double(n));

  for (int32_t bit = 0; bit < bits_; ++bit) {
    const int32_t offset = bins * (1 << bit);

#pragma omp simd aligned(idx1, idx2, d1, t : kAlign)
    // int64 operations are not vectorized in AVX2.
    for (int32_t i = 0; i < p; ++i) {
      const int64_t ts = int64_t(t[i]) + int64_t(offset);
      idx1[i] = MulAddPosMod(tf.a, int64_t(q) + int64_t(ts), tf.c, n);
      idx2[i] = MulAddPosMod(tf.a, int64_t(q) - int64_t(ts), tf.c, n);
      const int32_t k1 = MulMod(tf.b, ts, n);
      d1[i] = double(k1) / double(n) + PosModOne(delta * double(t[i]));
    }

#pragma omp simd aligned(idx1, idx2, s1, s2 : kAlign)
    // This doesn't seem vectorized either. Splitting into components helps.
    for (int32_t i = 0; i < p; ++i) {
      const Cplex x1 = x[idx1[i]] * bq;
      const Cplex x2 = std::conj(x[idx2[i]] * bq);
      s1[i] = 0.5 * (x1 + x2);
      s2[i] = RotateBackward(0.5 * (x1 - x2));
    }

    // CAUTION: Make sure the following is vectorized.
    for (int32_t i = 0; i < p; i += pack_card) {
      pack_t x(bs::aligned_load<pack_t>(d1_.data() + i));
      auto res = bs::sincos(2.0 * M_PI * x);
      bs::aligned_store(res.first, d2_.data() + i);
      bs::aligned_store(res.second, d3_.data() + i);
    }

#pragma omp simd aligned(s1, s2, d2, d3, wt : kAlign)
    for (int32_t i = 0; i < p; ++i) {
      const Cplex factor = Cplex(d3[i], d2[i]) * wt[i];
      s1[i] *= factor;
      s2[i] *= factor;
    }

    const int32_t folds = p / bins;
    for (int32_t i = 1; i < folds; ++i) {
      for (int32_t j = 0; j < bins; ++j) {
        const int32_t k = i * bins + j;
        s1[j] += s1[k];
        s2[j] += s2[k];
      }
    }

    // Do B-point FFT.
    plan_->Run(s1_, &(*out)[2 * bit + 1]);
    plan_->Run(s2_, &(*out)[2 * bit + 2]);

    // Combine the bin coefficients.
    for (int32_t bin = 0; bin < bins; ++bin) {
      const Cplex coef_r = (*out)[2 * bit + 1][bin];
      const Cplex coef_i = (*out)[2 * bit + 2][bin];
      (*out)[2 * bit + 1][bin] = coef_r + RotateForward(coef_i);
      (*out)[2 * bit + 2][bin] =
          std::conj(coef_r) + RotateForward(std::conj(coef_i));
    }
  }
}

BinInFreq::BinInFreq(const Window &win, int32_t bits)
    : win_(win), bits_(bits) {}

BinInFreq *BinInFreq::Create(int binner_type, const Window &win, int32_t bits) {
  if (binner_type == 0) {
    return new BinInFreqV0(win, bits);
  }
  if (binner_type == 1) {
    return new BinInFreqV1(win, bits);
  }
  LOG(FATAL) << "Unknown binner_type=" << binner_type;
  return nullptr;
}

BinInFreqV0::BinInFreqV0(const Window &win, int32_t bits)
    : BinInFreq(win, bits) {}

void BinInFreqV0::Run(const ModeMap &mm, const Transform &tf, int32_t q,
                      CplexMatrix *out) {
  CHECK_EQ(out->rows(), 1 + 2 * bits_);
  const int32_t bins = win_.bins();
  const int32_t n = win_.n();

  for (const auto &kv : mm) {
    const int32_t k0 = kv.first;
    const int32_t k1 = MulAddPosMod(tf.a, k0, tf.b, n); // 0 to n-1.
    const int32_t bin = int32_t(int64_t(k1) * int64_t(bins) / int64_t(n));
    const double xi =
        (double(bin) + 0.5) / double(bins) - double(k1) / double(n);
    const double wf = win_.SampleInFreq(xi);
    const Cplex base = kv.second * wf;
    for (int32_t u = 0; u < out->rows(); ++u) {
      const int32_t tau = GetTau(q, bins, u);
      const int32_t s = MulAddMod(tf.c, k0, int64_t(k1) * int64_t(tau), n);
      (*out)[u][bin] -= base * Sinusoid(double(s) / double(n));
    }
  }
}

BinInFreqV1::BinInFreqV1(const Window &win, int32_t bits)
    : BinInFreq(win, bits) {}

void BinInFreqV1::Run(const ModeMap &mm, const Transform &tf, int32_t q,
                      CplexMatrix *out) {
  CHECK_EQ(out->rows(), 1 + 2 * bits_);
  const int32_t bins = win_.bins();
  const int32_t n = win_.n();

  for (const auto &kv : mm) {
    const int32_t k0 = kv.first;
    const int32_t k1 = MulAddPosMod(tf.a, k0, tf.b, n);
    const int32_t bin = int32_t(int64_t(k1) * int64_t(bins) / int64_t(n));
    const double xi =
        (double(bin) + 0.5) / double(bins) - double(k1) / double(n);
    const double wf = win_.SampleInFreq(xi);

    const double freq =
        double(MulAddMod(q, k1, int64_t(tf.c) * int64_t(k0), n)) / double(n);
    const Cplex base = wf * kv.second * Sinusoid(freq);
    (*out)[0][bin] -= base;

    for (int32_t bit = 0; bit < bits_; ++bit) {
      const int32_t offset = bins * (1 << bit);
      const Cplex factor = Sinusoid(double(MulMod(offset, k1, n)) / double(n));
      (*out)[bit * 2 + 1][bin] -= base * factor;
      (*out)[bit * 2 + 2][bin] -= base * std::conj(factor); // Divide by factor.
    }
  }
}

} // namespace mps