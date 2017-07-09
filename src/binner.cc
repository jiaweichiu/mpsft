#include "binner.h"
#include "base.h"

namespace mps {

Binner::Binner(const BinnerOpt &opt) {
  CHECK_GT(opt.n, 0);
  CHECK_GT(opt.b, 0);
  CHECK_GT(opt.delta, 0);

  CHECK_EQ(1, opt.b % 2) << "Odd b expected";

  n_ = opt.n;
  b_ = opt.b;
  delta_ = opt.delta;

  width_ = 1.0 / (2.0 * b_);
  sqrt_c_delta_ = ::sqrt(-::log(delta_));
  sigma_f_ = 0.5 / (b_ * 2.0 * M_SQRT2 * sqrt_c_delta_);
  sigma_t_ = 1.0 / (2.0 * M_PI * sigma_f_);

  {
    // Decide p, the size of support.
    // p has to be sufficiently large.
    Real tmp = 2.0 * M_SQRT2 * sigma_t_ * sqrt_c_delta_ + 1;
    Int factor = Round(tmp / Real(b_));
    if ((factor % 2) == 0) {
      ++factor;
    }
    // p has to be odd and a multiple of b.
    p_ = factor * b_;
  }
  LOG(INFO) << b_ << " " << p_;
}

} // namespace mps