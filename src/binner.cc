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

// delta = -0.5/bins.
// bq_factor = Sinusoid(b*q/N).
// win_factor = wt * 0.5.
// inline CplexPair BinInTimeHelper(Int n, double delta, const Transform &tf,
// Int
// t,
//                                  Int q, Int s, Cplex bq_factor, double
//                                  win_factor,
//                                  const vector<Cplex> &x) {
//   const Int i = PosMod(Long(tf.a) * Long(q + s + t) + Long(tf.c), n);
//   const Int j = PosMod(Long(tf.a) * Long(q - s - t) + Long(tf.c), n);
//   const Int u = Mod(Long(tf.b) * Long(t + s), n);
//   const Cplex sinusoid = Sinusoid(double(u) / double(n) - delta);
//   const Cplex v1 = x[i] * bq_factor;
//   const Cplex v2 = std::conj(x[j] * bq_factor);

//   return std::make_pair(win_factor * sinusoid * (v1 + v2),
//                         RotateBackward(win_factor * sinusoid * (v1 - v2)));
// }

Binner::Binner(const Window &win, Int bits) : win_(win), bits_(bits) {
  const Int bins = win.bins();
  plan_.reset(new FFTPlan(bins, FFTW_FORWARD));
  scratch_.reset(new CplexArray(bins));
}

BinnerSimple::BinnerSimple(const Window &win, Int bits) : Binner(win, bits) {}

void BinnerSimple::BinInTime(const CplexArray &x, const Transform &tf, Int q,
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
    scratch_->Clear();
    for (Int i = 0; i < p; ++i) {
      const Int t = i <= p2 ? i : i - p;
      const double wt = win_.wt(i <= p2 ? i : p - i);
      const Int j = PosMod(Long(tf.a) * Long(t + tau) + Long(tf.c), n);
      const Int k = Mod(Long(tf.b) * Long(t + tau), n);
      // Do mods for better accuracy. Note fmod can be negative, but it is ok.
      const double angle = (2.0 * M_PI) * (double(k) / double(n) +
                                           std::fmod(delta * double(t), 1.0));
      (*scratch_)[i % bins] += (x[j] * Sinusoid(angle)) * wt;
    }
    // Do B-point FFT.
    CplexArray &v = (*out)[u];
    DCHECK_EQ(bins, v.size());
    plan_->Run(*scratch_, &v);
  }
}

void BinnerSimple::BinInFreq(const ModeMap &mm, const Transform &tf, Int q,
                             CplexMatrix *out) {
  CHECK_EQ(out->rows(), 1 + 2 * bits_);
  const Int bins = win_.bins();
  const Int n = win_.n();

  for (const auto &kv : mm) {
    const Int k = kv.first;
    const Int l = PosMod(Long(tf.a) * Long(k) + Long(tf.b), n); // 0 to n-1.
    const Int bin = Int(Long(l) * Long(bins) / Long(n));
    const double xi =
        (double(bin) + 0.5) / double(bins) - double(l) / double(n);
    const double wf = win_.SampleInFreq(xi);
    const Cplex base = kv.second * wf;
    for (Int u = 0; u < out->rows(); ++u) {
      const Int tau = GetTau(q, bins, u);
      const Int s = Mod(Long(tf.c) * Long(k) + Long(l) * Long(tau), n);
      const double angle = (2.0 * M_PI) * (double(s) / double(n));
      (*out)[u][bin] -= base * Sinusoid(angle);
    }
  }
}

BinnerFast::BinnerFast(const Window &win, Int bits) : Binner(win, bits) {
  scratch2_.reset(new CplexArray(win.bins()));
}

void BinnerFast::BinInTime(const CplexArray &x, const Transform &tf, Int q,
                           CplexMatrix *out) {
  CHECK_EQ(out->rows(), 1 + 2 * bits_);
  const Int bins = win_.bins();
  const Int n = win_.n();

  const double delta = -0.5 / double(bins);
  DCHECK_EQ(n, x.size());

  const Int p = win_.p();
  const Int p2 = (p - 1) / 2;

  // Take care of tau=q first.
  scratch_->Clear();
  for (Int i = 0; i < p; ++i) {
    const Int t = i <= p2 ? i : i - p;
    const double wt = win_.wt(i <= p2 ? i : p - i);
    const Int j = PosMod(Long(tf.a) * Long(t + q) + Long(tf.c), n);
    const Int k = Mod(Long(tf.b) * Long(t + q), n);
    const double angle = (2.0 * M_PI) * (double(k) / double(n) +
                                         std::fmod(delta * double(t), 1.0));
    (*scratch_)[i % bins] += (x[j] * Sinusoid(angle)) * wt;
  }
  plan_->Run(*scratch_, &(*out)[0]);

  // Take care of offsets from q.
  const Cplex bq =
      Sinusoid((2.0 * M_PI) * double(Mod(Long(tf.b) * Long(q), n)) / double(n));

  for (Int bit = 0; bit < bits_; ++bit) {
    const Int offset = bins * (1 << bit);
    scratch_->Clear();
    scratch2_->Clear();
    for (Int i = 0; i < p; ++i) {
      const Int t = i <= p2 ? i : i - p;
      const double wt = win_.wt(i <= p2 ? i : p - i);

      const Long ts = Long(t) + Long(offset);
      const Int j1 = PosMod(Long(tf.a) * (Long(q) + ts) + Long(tf.c), n);
      const Int j2 = PosMod(Long(tf.a) * (Long(q) - ts) + Long(tf.c), n);
      const Cplex x1 = x[j1];
      const Cplex x2 = x[j2];

      const Int k = Mod(Long(tf.b) * ts, n);
      // Do mods for better accuracy. Note fmod can be negative, but it is ok.
      const double angle = (2.0 * M_PI) * (double(k) / double(n) +
                                           std::fmod(delta * double(t), 1.0));
      // Only one Sinusoid per iteration.
      const Cplex factor = Sinusoid(angle);

      const Cplex z1 = x1 * bq * factor * wt;
      const Cplex z2 = std::conj(x2) * std::conj(bq) * factor * wt;

      (*scratch_)[i % bins] += 0.5 * (z1 + z2);
      (*scratch2_)[i % bins] += RotateBackward(0.5 * (z1 - z2));
    }
    // Do B-point FFT.
    plan_->Run(*scratch_, &(*out)[2 * bit + 1]);
    plan_->Run(*scratch2_, &(*out)[2 * bit + 2]);

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

void BinnerFast::BinInFreq(const ModeMap &mm, const Transform &tf, Int q,
                           CplexMatrix *out) {
  CHECK_EQ(out->rows(), 1 + 2 * bits_);
  const Int bins = win_.bins();
  const Int n = win_.n();

  for (const auto &kv : mm) {
    const Int k0 = kv.first;
    const Int k1 = PosMod(Long(tf.a) * Long(k0) + Long(tf.b), n); // 0 to n-1.
    const Int bin = Int(Long(k1) * Long(bins) / Long(n));
    const double xi =
        (double(bin) + 0.5) / double(bins) - double(k1) / double(n);
    const double wf = win_.SampleInFreq(xi);

    const double angle =
        (2.0 * M_PI) *
        double(Mod(Long(q) * Long(k1) + Long(tf.c) * Long(k0), n)) / double(n);
    const Cplex base = wf * kv.second * Sinusoid(angle);
    (*out)[0][bin] -= base;

    for (Int bit = 0; bit < bits_; ++bit) {
      const Int offset = bins * (1 << bit);
      const Cplex factor = Sinusoid(
          (2.0 * M_PI) * double(Mod(Long(offset) * Long(k1), n)) / double(n));
      (*out)[bit * 2 + 1][bin] -= base * factor;
      (*out)[bit * 2 + 2][bin] -= base * std::conj(factor); // Divide by factor.
    }
  }
}

Binner *Binner::Create(int binner_type, const Window &win, Int bits) {
  switch (binner_type) {
  case kBinnerSimple: {
    return new BinnerSimple(win, bits);
  }
  case kBinnerFast: {
    return new BinnerFast(win, bits);
  }
  default: { LOG(FATAL) << "Unknown binner_type: " << binner_type; }
  }
  return nullptr;
}

} // namespace mps