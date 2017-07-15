#include "binner.h"
#include "base.h"
#include "window.h"

namespace mps {

TauSet::TauSet(Int q, Int bins, Int bits) : q_(q), bins_(bins), bits_(bits) {}

Int TauSet::tau(Int idx) const {
  if (idx == 0) {
    return q_;
  }
  if (idx & 1) {
    return q_ + bins_ * (1 << ((idx - 1) / 2));
  }
  return q_ - bins_ * (1 << (idx / 2 - 1));
}

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

void BinInTime(const Window &win, const Transform &tf, const TauSet &taus,
               const CplexArray &x, FFTPlan *plan, CplexMatrix *out,
               CplexArray *scratch) {
  // const Int q = taus.q;
  DCHECK(plan);
  DCHECK(scratch);
  DCHECK_EQ(out->rows(), taus.size());
  DCHECK_EQ(plan->sign(), FFTW_FORWARD);

  const Int bins = win.bins();
  DCHECK_EQ(bins, scratch->size());
  DCHECK_EQ(bins, plan->n());

  const double delta = -0.5 / double(bins);
  // const Cplex bq_factor = Sinusoid(Mod(Long(tf.b) * Long(q), n));

  const Int n = win.n();
  DCHECK_EQ(n, x.size());

  const Int p = win.p();
  const Int p2 = (p - 1) / 2;

  for (Int u = 0; u < taus.size(); ++u) {
    const Int tau = taus.tau(u);
    scratch->Clear();
    for (Int i = 0; i < p; ++i) {
      const Int t = i <= p2 ? i : i - p;
      const double wt = win.wt(i <= p2 ? i : p - i);
      const Int j = PosMod(Long(tf.a) * Long(t + tau) + Long(tf.c), n);
      const Int k = PosMod(Long(tf.b) * Long(t + tau), n);
      // Do mods for better accuracy. Note fmod can be negative, but it is ok.
      const double angle = (2.0 * M_PI) * (double(k) / double(n) +
                                           std::fmod(delta * double(t), 1.0));
      (*scratch)[i % bins] += (x[j] * Sinusoid(angle)) * wt;
    }
    // Do B-point FFT.
    CplexArray &v = (*out)[u];
    DCHECK_EQ(bins, v.size());
    plan->Run(*scratch, &v);
  }
}

void BinInFreq(const Window &win, const Transform &tf, const TauSet &taus,
               const ModeMap &mm, CplexMatrix *out) {
  const Int bins = win.bins();
  const Int n = win.n();

  for (const auto &kv : mm) {
    const Int k = kv.first;
    const Int l = PosMod(Long(tf.a) * Long(k) + Long(tf.b), n); // 0 to n-1.
    const Int bin = Int(Long(l) * Long(bins) / Long(n));
    const double xi =
        (double(bin) + 0.5) / double(bins) - double(l) / double(n);
    const double wf = win.SampleInFreq(xi);
    for (Int u = 0; u < taus.size(); ++u) {
      const Int tau = taus.tau(u);
      const Int s = Mod(Long(tf.c) * Long(k) + Long(l) * Long(tau), n);
      const double angle = (2.0 * M_PI) * (double(s) / double(n));
      (*out)[u][bin] -= (kv.second * Sinusoid(angle)) * wf;
    }
  }
}

} // namespace mps