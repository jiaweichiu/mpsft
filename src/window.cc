/*
 * Copyright (c) 2017 Jiawei Chiu
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301 USA
 *
 */
#include <boost/math/special_functions/erf.hpp>
#include <boost/math/special_functions/sinc.hpp>
#include <glog/logging.h>

#include "base.h"
#include "window.h"

namespace mps {

Window::Window(Int n, Int bins, double delta)
    : n_(n), bins_(bins), delta_(delta) {
  CHECK_GT(n_, 0);
  CHECK_GT(bins_, 0);
  CHECK_GT(delta_, 0);
  CHECK_EQ(1, bins_ % 2) << "Odd b expected";

  width_ = 1.0 / (2.0 * bins_);
  const double sqrt_c_delta = ::sqrt(-::log(delta_));
  sigma_f_ = 0.5 / (bins_ * 2.0 * M_SQRT2 * sqrt_c_delta);
  sigma_t_ = 1.0 / ((2.0 * M_PI) * sigma_f_);

  {
    // Decide p, the size of support.
    // p has to be sufficiently large.
    double tmp = 2.0 * M_SQRT2 * sigma_t_ * sqrt_c_delta + 1;
    Int factor = Int(std::ceil(tmp / double(bins_)));
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

double Window::SampleInTime(Int i) const {
  const double t = double(i);
  const double u = t * M_PI * sigma_f_;
  return width_ * boost::math::sinc_pi<double>(t * M_PI * width_) *
         ::exp(-2.0 * u * u);
}

double Window::SampleInFreq(double xi) const {
  const double c = 0.5 * width_;
  const double d = M_SQRT2 * M_PI * sigma_t_;
  return 0.5 * (boost::math::erf<double>((xi + c) * d) -
                boost::math::erf<double>((xi - c) * d));
}

} // namespace mps