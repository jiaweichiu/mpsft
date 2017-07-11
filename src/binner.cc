#include "binner.h"
#include "base.h"
#include "window.h"

namespace mps {

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
               const CplexArray &x, FFTPlan *plan, CplexArray *out1,
               CplexArray *out2) {
  // const Int q = taus.q;
  const Int bins = win.bins();
  DCHECK_EQ(bins, out1->Size());
  DCHECK_EQ(bins, out2->Size());
  DCHECK_EQ(bins, plan->n());

  const Real delta = -0.5 / Real(bins);
  // const Cplex bq_factor = Sinusoid(Mod(Long(tf.b) * Long(q), n));

  const Int n = win.n();
  DCHECK_EQ(n, x.Size());

  const Int p = win.p();
  const Int p2 = (p - 1) / 2;

  const Int tau = taus.q; // CHANGE

  out1->Clear();

  for (Int i = 0; i < p; ++i) {
    const Int t = i <= p2 ? i : i - p;
    const Real wt = win.wt(i <= p2 ? i : p - i);
    // LOG(INFO) << "t=" << t << " wt=" << wt;
    const Int j = PosMod(Long(tf.a) * Long(t + tau) + Long(tf.c), n);
    const Int k = PosMod(Long(tf.b) * Long(t + tau), n);
    const Real angle = (2.0 * M_PI) * (Real(k) / Real(n) - delta);
    (*out1)[i % bins] += (x[j] * Sinusoid(angle)) * wt;
  }
  // Do B-point FFT.
  plan->Run(*out1, out2);
}

void BinInFreq(const Window &win, const Transform &tf, const TauSet &taus,
               const CplexArray &coef, const vector<Int> &loc,
               CplexArray *out) {
  const Int tau = taus.q; // CHANGE.

  DCHECK_EQ(coef.Size(), loc.size());
  const Int bins = win.bins();
  const Int n = win.n();
  for (Int i = 0; i < coef.Size(); ++i) {
    const Int k = loc[i];
    const Int l = PosMod(Long(tf.a) * Long(k) + Long(tf.b), n); // 0 to n-1.
    const Int bin = Int(Long(l) * Long(bins) / Long(n));
    const Real xi = (Real(bin) + 0.5) / Real(bins) - Real(l) / Real(n);
    const Real wf = win.SampleInFreq(xi);
    const Int s = Mod(Long(tf.c) * Long(k) + Long(l) * Long(tau), n);
    const Real angle = (2.0 * M_PI) * (Real(s) / Real(n));
    (*out)[i] -= (coef[i] * Sinusoid(angle)) * wf;
  }
}

} // namespace mps