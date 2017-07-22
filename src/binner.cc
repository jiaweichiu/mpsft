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
#include "integer.h"
#include "sincos.h"
#include "window.h"

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
  // if (binner_type == 1) {
  //   return new BinInTimeV1(win, bits);
  // }
  // if (binner_type == 2) {
  //   return new BinInTimeV2(win, bits);
  // }
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
      const double wt = win_.wt(i);
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

// BinInTimeV1::BinInTimeV1(const Window &win, int32_t bits)
//     : BinInTime(win, bits), scratch_(win.bins()), scratch2_(win.bins()) {}

// inline void BinInTimeV1Helper(int32_t i, int32_t t, double wt, int32_t
// offset, int32_t q, int32_t n,
//                               int32_t bins, Cplex bq, double delta,
//                               const Transform &tf, const CplexArray &x,
//                               CplexArray *scratch, CplexArray *scratch2) {
//   const int32_t ts = t + offset;
//   const int32_t j1 = PosMod(tf.a * (q + ts) + tf.c, n);
//   const int32_t j2 = PosMod(tf.a * (q - ts) + tf.c, n);
//   const Cplex x1 = x[j1];
//   const Cplex x2 = x[j2];

//   const int32_t k = Mod(tf.b * ts, n);
//   // Do mods for better accuracy. Note fmod can be negative, but it is ok.
//   const double freq = double(k) / double(n) + std::fmod(delta * double(t),
//   1.0);
//   // Only one Sinusoid per iteration.
//   const Cplex factor = Sinusoid(freq);

//   const Cplex z1 = x1 * bq * factor * wt;
//   const Cplex z2 = std::conj(x2) * std::conj(bq) * factor * wt;

//   (*scratch)[i % bins] += 0.5 * (z1 + z2);
//   (*scratch2)[i % bins] += RotateBackward(0.5 * (z1 - z2));
// }

// void BinInTimeV1::Run(const CplexArray &x, const Transform &tf, int32_t q,
//                       CplexMatrix *out) {
//   CHECK_EQ(out->rows(), 1 + 2 * bits_);
//   const int32_t bins = win_.bins();
//   const int32_t n = win_.n();

//   const double delta = -0.5 / double(bins);
//   DCHECK_EQ(n, x.size());

//   const int32_t p = win_.p();
//   const int32_t p2 = (p - 1) / 2;

//   // Take care of tau=q first.
//   scratch_.Clear();
//   for (int32_t i = 0; i < p; ++i) {
//     const int32_t t = i <= p2 ? i : i - p;
//     const double wt = win_.wt(i <= p2 ? i : p - i);
//     const int32_t j = PosMod(tf.a * (t + q) + tf.c, n);
//     const int32_t k = Mod(tf.b * (t + q), n);
//     const double freq =
//         double(k) / double(n) + std::fmod(delta * double(t), 1.0);
//     scratch_[i % bins] += (x[j] * Sinusoid(freq)) * wt;
//   }
//   plan_->Run(scratch_, &(*out)[0]);

//   // Take care of offsets from q.
//   const Cplex bq = Sinusoid(double(Mod(tf.b * q, n)) / double(n));

//   for (int32_t bit = 0; bit < bits_; ++bit) {
//     const int32_t offset = bins * (1 << bit);
//     scratch_.Clear();
//     scratch2_.Clear();
//     for (int32_t i = 0; i <= p2; ++i) {
//       const int32_t t = i;
//       const double wt = win_.wt(i);
//       BinInTimeV1Helper(i, t, wt, offset, q, n, bins, bq, delta, tf, x,
//                         &scratch_, &scratch2_);
//     }
//     for (int32_t i = p2 + 1; i < p; ++i) {
//       const int32_t t = i - p;
//       const double wt = win_.wt(p - i);
//       BinInTimeV1Helper(i, t, wt, offset, q, n, bins, bq, delta, tf, x,
//                         &scratch_, &scratch2_);
//     }

//     // Do B-point FFT.
//     plan_->Run(scratch_, &(*out)[2 * bit + 1]);
//     plan_->Run(scratch2_, &(*out)[2 * bit + 2]);

//     // Combine the bin coefficients.
//     for (int32_t bin = 0; bin < bins; ++bin) {
//       const Cplex coef_r = (*out)[2 * bit + 1][bin];
//       const Cplex coef_i = (*out)[2 * bit + 2][bin];
//       (*out)[2 * bit + 1][bin] = coef_r + RotateForward(coef_i);
//       (*out)[2 * bit + 2][bin] =
//           std::conj(coef_r) + RotateForward(std::conj(coef_i));
//     }
//   }
// }

// BinInTimeV2::BinInTimeV2(const Window &win, int32_t bits)
//     : BinInTime(win, bits), scratch_(win.p()), scratch2_(win.p()),
//       idx_(win.p()), idx2_(win.p()) {}

// inline void BinInTimeV2Helper1(int32_t i, int32_t t, const Transform &tf,
// int32_t offset,
//                                int32_t q, int32_t n, int32_t *idx, int32_t
//                                *idx2) {
//   const int32_t ts = t + offset;
//   idx[i] = PosMod(tf.a * (q + ts) + tf.c, n);
//   idx2[i] = PosMod(tf.a * (q - ts) + tf.c, n);
// }

// inline void BinInTimeV2Helper2(int32_t i, Cplex bq, const CplexArray &x,
// int32_t *idx,
//                                int32_t *idx2, Cplex *out, Cplex *out2) {
//   const Cplex x1 = x[idx[i]] * bq;
//   const Cplex x2 = std::conj(x[idx2[i]] * bq);
//   out[i] = 0.5 * (x1 + x2);
//   out2[i] = RotateBackward(0.5 * (x1 - x2));
// }

// inline void BinInTimeV2Helper3(int32_t i, int32_t t, const Transform &tf,
// int32_t offset,
//                                double wt, int32_t n, double delta, Cplex
//                                *out,
//                                Cplex *out2) {
//   const int32_t ts = t + offset;
//   const int32_t k = Mod(tf.b * ts, n);
//   // Do mods for better accuracy. Note fmod can be negative, but it is ok.
//   const double freq = double(k) / double(n) + std::fmod(delta * double(t),
//   1.0);
//   const Cplex factor = Sinusoid(freq) * wt;
//   out[i] *= factor;
//   out2[i] *= factor;
// }

// void BinInTimeV2::Run(const CplexArray &x, const Transform &tf, int32_t q,
//                       CplexMatrix *out) {
//   CHECK_EQ(out->rows(), 1 + 2 * bits_);
//   const int32_t bins = win_.bins();
//   const int32_t n = win_.n();

//   const double delta = -0.5 / double(bins);
//   DCHECK_EQ(n, x.size());

//   const int32_t p = win_.p();
//   const int32_t p2 = (p - 1) / 2;

//   // Take care of tau=q first. Treat this case simply.
//   scratch_.Clear();
//   for (int32_t i = 0; i < p; ++i) {
//     const int32_t t = i <= p2 ? i : i - p;
//     const double wt = win_.wt(i <= p2 ? i : p - i);
//     const int32_t j = PosMod(tf.a * (t + q) + tf.c, n);
//     const int32_t k = Mod(tf.b * (t + q), n);
//     const double freq =
//         double(k) / double(n) + std::fmod(delta * double(t), 1.0);
//     scratch_[i % bins] += (x[j] * Sinusoid(freq)) * wt;
//   }
//   plan_->Run(scratch_, &(*out)[0]);

//   // Take care of offsets from q.
//   const Cplex bq = Sinusoid(double(Mod(tf.b * q, n)) / double(n));

//   for (int32_t bit = 0; bit < bits_; ++bit) {
//     const int32_t offset = bins * (1 << bit);
//     Cplex *__restrict__ s1 = scratch_.data();
//     Cplex *__restrict__ s2 = scratch2_.data();
//     int32_t *__restrict__ idx1 = idx_.data();
//     int32_t *__restrict__ idx2 = idx2_.data();

// #pragma omp simd aligned(idx1, idx2 : kAlign)
//     for (int32_t i = 0; i <= p2; ++i) {
//       BinInTimeV2Helper1(i, i, tf, offset, q, n, idx1, idx2);
//     }

// #pragma omp simd aligned(idx1, idx2 : kAlign)
//     for (int32_t i = p2 + 1; i < p; ++i) {
//       BinInTimeV2Helper1(i, i - p, tf, offset, q, n, idx1, idx2);
//     }

// #pragma omp simd aligned(idx1, idx2, s1, s2 : kAlign)
//     for (int32_t i = 0; i < p; ++i) {
//       BinInTimeV2Helper2(i, bq, x, idx1, idx2, s1, s2);
//     }

// #pragma omp simd aligned(s1, s2 : kAlign)
//     for (int32_t i = 0; i <= p2; ++i) {
//       BinInTimeV2Helper3(i, i, tf, offset, win_.wt(i), n, delta, s1, s2);
//     }

// #pragma omp simd aligned(s1, s2 : kAlign)
//     for (int32_t i = p2 + 1; i < p; ++i) {
//       BinInTimeV2Helper3(i, i - p, tf, offset, win_.wt(p - i), n, delta, s1,
//                          s2);
//     }

//     const int32_t folds = p / bins;
//     for (int32_t i = 1; i < folds; ++i) {
//       for (int32_t j = 0; j < bins; ++j) {
//         const int32_t k = i * bins + j;
//         scratch_[j] += scratch_[k];
//         scratch2_[j] += scratch2_[k];
//       }
//     }

//     // Do B-point FFT.
//     plan_->Run(scratch_, &(*out)[2 * bit + 1]);
//     plan_->Run(scratch2_, &(*out)[2 * bit + 2]);

//     // Combine the bin coefficients.
//     for (int32_t bin = 0; bin < bins; ++bin) {
//       const Cplex coef_r = (*out)[2 * bit + 1][bin];
//       const Cplex coef_i = (*out)[2 * bit + 2][bin];
//       (*out)[2 * bit + 1][bin] = coef_r + RotateForward(coef_i);
//       (*out)[2 * bit + 2][bin] =
//           std::conj(coef_r) + RotateForward(std::conj(coef_i));
//     }
//   }
// }

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