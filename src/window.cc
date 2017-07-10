#include <boost/math/special_functions/erf.hpp>
#include <boost/math/special_functions/sinc.hpp>
#include <glog/logging.h>

#include "base.h"
#include "window.h"

namespace mps {

Window::Window(Int n, Int bins, Real delta)
    : n_(n), bins_(bins), delta_(delta) {
  CHECK_GT(n_, 0);
  CHECK_GT(bins_, 0);
  CHECK_GT(delta_, 0);
  CHECK_EQ(1, bins_ % 2) << "Odd b expected";

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
    wt_[i] = SampleInTime(i);
  }
}

Real Window::SampleInTime(Int i) const {
  const Real t = Real(i);
  const Real u = t * M_PI * sigma_f_;
  return width_ * boost::math::sinc_pi<Real>(t * M_PI * width_) *
         ::exp(-2.0 * u * u);
}

Real Window::SampleInFreq(Real xi) const {
  const Real c = 0.5 * width_;
  // const Real d = 1.0 / (M_SQRT2 * sigma_f_);
  const Real d = M_SQRT2 * M_PI * sigma_t_;
  return 0.5 * (boost::math::erf<Real>((xi + c) * d) -
                boost::math::erf<Real>((xi - c) * d));
}

// Real Window::Energy() const {
//   const Int p2 = (p_ - 1) / 2;
//   Real out = Square(wt_[0]);
//   for (Int i = 1; i <= p2; ++i) {
//     out += 2 * Square(wt_[i]);
//   }
//   return out;
// }

} // namespace mps