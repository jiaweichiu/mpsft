#include "binner.h"
#include "base.h"
#include "window.h"

namespace mps {

Int TauSet::value(Int idx) const {
  if (idx == 0) {
    return q;
  }
  if (idx & 1) {
    return q + list_s[(idx - 1) / 2];
  }
  return q - list_s[idx / 2 - 1];
}

// delta = -0.5/bins.
// bq_factor = Sinusoid(b*q/N).
// win_factor = wt * 0.5.
// inline CplexPair BinInTimeHelper(Int n, Real delta, const Transform &tf, Int
// t,
//                                  Int q, Int s, Cplex bq_factor, Real
//                                  win_factor,
//                                  const vector<Cplex> &x) {
//   const Int i = PosMod(Long(tf.a) * Long(q + s + t) + Long(tf.c), n);
//   const Int j = PosMod(Long(tf.a) * Long(q - s - t) + Long(tf.c), n);
//   const Int u = Mod(Long(tf.b) * Long(t + s), n);
//   const Cplex sinusoid = Sinusoid(Real(u) / Real(n) - delta);
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

  const Real delta = -0.5 / Real(bins);
  // const Cplex bq_factor = Sinusoid(Mod(Long(tf.b) * Long(q), n));

  const Int n = win.n();
  DCHECK_EQ(n, x.size());

  const Int p = win.p();
  const Int p2 = (p - 1) / 2;

  for (Int u = 0; u < taus.size(); ++u) {
    const Int tau = taus.value(u);
    scratch->Clear();
    for (Int i = 0; i < p; ++i) {
      const Int t = i <= p2 ? i : i - p;
      const Real wt = win.wt(i <= p2 ? i : p - i);
      const Int j = PosMod(Long(tf.a) * Long(t + tau) + Long(tf.c), n);
      const Int k = PosMod(Long(tf.b) * Long(t + tau), n);
      // Do mods for better accuracy. Note fmod can be negative, but it is ok.
      const Real angle =
          (2.0 * M_PI) * (Real(k) / Real(n) + std::fmod(delta * Real(t), 1.0));
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

  for (const auto& kv : mm) {
    const Int k = kv.first;
    const Int l = PosMod(Long(tf.a) * Long(k) + Long(tf.b), n); // 0 to n-1.
    const Int bin = Int(Long(l) * Long(bins) / Long(n));
    const Real xi = (Real(bin) + 0.5) / Real(bins) - Real(l) / Real(n);
    const Real wf = win.SampleInFreq(xi);
    for (Int u = 0; u < taus.size(); ++u) {
      const Int tau = taus.value(u);
      const Int s = Mod(Long(tf.c) * Long(k) + Long(l) * Long(tau), n);
      const Real angle = (2.0 * M_PI) * (Real(s) / Real(n));
      (*out)[u][bin] -= (kv.second * Sinusoid(angle)) * wf;
    }
  }
}

} // namespace mps