#include <boost/math/special_functions/erf.hpp>
#include <boost/math/special_functions/sinc.hpp>

#include "base.h"
#include "binner.h"

namespace mps {

Window::Window(const BinnerOpt& opt) {
  CHECK_GT(opt.n, 0);
  CHECK_GT(opt.bins, 0);
  CHECK_GT(opt.delta, 0);

  CHECK_EQ(1, opt.bins % 2) << "Odd b expected";

  n_ = opt.n;
  bins_ = opt.bins;
  delta_ = opt.delta;

  width_ = 1.0 / (2.0 * bins_);
  sqrt_c_delta_ = ::sqrt(-::log(delta_));
  sigma_f_ = 0.5 / (bins_ * 2.0 * M_SQRT2 * sqrt_c_delta_);
  sigma_t_ = 1.0 / (2.0 * M_PI * sigma_f_);

  {
    // Decide p, the size of support.
    // p has to be sufficiently large.
    Real tmp = 2.0 * M_SQRT2 * sigma_t_ * sqrt_c_delta_ + 1;
    Int factor = Round(tmp / Real(bins_));
    if ((factor % 2) == 0) {
      ++factor;
    }
    // p has to be odd and a multiple of b.
    p_ = factor * bins_;
  }
  CHECK_EQ(1, p_ % 2) << "Odd p expected";

  const Int p2 = (p_ - 1) / 2;
  wt_.resize(p2 + 1);
  wt_[0] = width_;
  for (Int i = 1; i <= p2; ++i) {
    const Real t = Real(i);
    const Real u = t * M_PI * sigma_f_;
    const double v = width_ * boost::math::sinc_pi<Real>(t * M_PI * width_) *
                     ::exp(-2.0 * u * u);
    wt_[i] = v;
  }
}

// delta = -0.5/bins.
// bq_factor = Sinusoid(b*q/N).
// win_factor = wt * 0.5.
inline CplexPair BinInTimeHelper(Int n, Real delta, const Transform& tf, Int t,
                                 Int q, Int s, Cplex bq_factor, Real win_factor,
                                 const vector<Cplex>& x) {
  const Int i = PosMod(Long(tf.a) * Long(q + s + t) + Long(tf.c), n);
  const Int j = PosMod(Long(tf.a) * Long(q - s - t) + Long(tf.c), n);
  const Int u = Mod(Long(tf.b) * Long(t + s), n);
  const Cplex sinusoid = Sinusoid(Real(u) / Real(n) - delta);
  const Cplex v1 = x[i] * bq_factor;
  const Cplex v2 = std::conj(x[j] * bq_factor);

  return std::make_pair(win_factor * sinusoid * (v1 + v2),
                        RotateBackward(win_factor * sinusoid * (v1 - v2)));
}

void BinInTime(const Window& win, const Transform& tf, const TauSet& taus,
               const vector<Cplex>& x) {}

}  // namespace mps