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
#include <cmath>
#include <glog/logging.h>

#include "base.h"
#include "sincos.h"
#include "window.h"

namespace mps {

Window::Window(int32_t n, int32_t bins, double delta)
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
    int32_t factor = int32_t(std::ceil(tmp / double(bins_)));
    if ((factor % 2) == 0) {
      ++factor;
    }
    // p has to be odd and a multiple of b.
    p_ = factor * bins_;
  }
  CHECK_EQ(1, p_ % 2) << "Odd p expected";

  const int32_t p2 = (p_ - 1) / 2;

  t_.resize(p_);
  t_[0] = 0;
  for (int32_t i = 1; i <= p2; ++i) {
    t_[i] = i;
  }
  for (int32_t i = p2 + 1; i < p_; ++i) {
    t_[i] = i - p_;
  }

  wt_.resize(p_);
  wt_[0] = width_;
  for (int32_t i = 1; i < p_; ++i) {
    wt_[i] = SampleInTime(t_[i]);
  }
}

double Window::SampleInTime(int32_t i) const {
  const double t = double(i);
  const double u = t * M_PI * sigma_f_;
  return width_ * SincPi(t * M_PI * width_) * std::exp(-2.0 * u * u);
}

double Window::SampleInFreq(double xi) const {
  const double c = 0.5 * width_;
  const double d = M_SQRT2 * M_PI * sigma_t_;
  return 0.5 * (std::erf((xi + c) * d) - std::erf((xi - c) * d));
}

} // namespace mps