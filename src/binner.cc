#include "binner.h"
#include "base.h"

namespace mps {

// delta = -0.5/bins.
// bq_factor = Sinusoid(b*q/N).
// win_factor = wt * 0.5.
inline CplexPair BinInTimeHelper(Int n, Real delta, const Transform &tf, Int t,
                                 Int q, Int s, Cplex bq_factor, Real win_factor,
                                 const vector<Cplex> &x) {
  const Int i = PosMod(Long(tf.a) * Long(q + s + t) + Long(tf.c), n);
  const Int j = PosMod(Long(tf.a) * Long(q - s - t) + Long(tf.c), n);
  const Int u = Mod(Long(tf.b) * Long(t + s), n);
  const Cplex sinusoid = Sinusoid(Real(u) / Real(n) - delta);
  const Cplex v1 = x[i] * bq_factor;
  const Cplex v2 = std::conj(x[j] * bq_factor);

  return std::make_pair(win_factor * sinusoid * (v1 + v2),
                        RotateBackward(win_factor * sinusoid * (v1 - v2)));
}

void BinInTime(const Window &win, const Transform &tf, const TauSet &taus,
               const vector<Cplex> &x) {
  const Int q = taus.q;
  const Int bins = win.bins();
  const Real delta = -0.5 / Real(bins);
  // const Cplex bq_factor = Sinusoid(Mod(Long(tf.b) * Long(q), n));
}

} // namespace mps