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

void Binner::BinInTime(const CplexArray &x, Int q, CplexMatrix *out) {
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
      const Int j = PosMod(Long(tf_.a) * Long(t + tau) + Long(tf_.c), n);
      const Int k = PosMod(Long(tf_.b) * Long(t + tau), n);
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

void Binner::BinInFreq(const ModeMap &mm, Int q, CplexMatrix *out) {
  CHECK_EQ(out->rows(), 1 + 2 * bits_);
  const Int bins = win_.bins();
  const Int n = win_.n();

  for (const auto &kv : mm) {
    const Int k = kv.first;
    const Int l = PosMod(Long(tf_.a) * Long(k) + Long(tf_.b), n); // 0 to n-1.
    const Int bin = Int(Long(l) * Long(bins) / Long(n));
    const double xi =
        (double(bin) + 0.5) / double(bins) - double(l) / double(n);
    const double wf = win_.SampleInFreq(xi);
    for (Int u = 0; u < out->rows(); ++u) {
      const Int tau = GetTau(q, bins, u);
      const Int s = Mod(Long(tf_.c) * Long(k) + Long(l) * Long(tau), n);
      const double angle = (2.0 * M_PI) * (double(s) / double(n));
      (*out)[u][bin] -= (kv.second * Sinusoid(angle)) * wf;
    }
  }
}

Binner::Binner(const Window &win, const Transform &tf, Int bits)
    : win_(win), tf_(tf), bits_(bits) {
  const Int bins = win.bins();
  plan_.reset(new FFTPlan(bins, FFTW_FORWARD));
  scratch_.reset(new CplexArray(bins));
}

} // namespace mps