#include <cmath>
#include <glog/logging.h>

#include "binner.h"

namespace mps {

Binner::Binner(const BinnerOpt &opt) {
  CHECK_GT(opt.n, 0);
  CHECK_GT(opt.b, 0);
  CHECK_GT(opt.delta, 0);

  n_ = opt.n;
  b_ = opt.b;
  delta_ = opt.delta;

  width_ = 1.0 / (2.0 * b_);
  sqrt_c_delta_ = ::sqrt(-::log(delta_));
  sigma_f_ = 0.5 / (b_ * 2 * M_SQRT2 * sqrt_c_delta_);
  sigma_t_ = 1.0 / (2 * M_PI * sigma_f_);
}

}  // namespace mps