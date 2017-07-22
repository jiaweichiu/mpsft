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
#pragma once

#include "base.h"

namespace mps {

class Window {
public:
  Window(int32_t n, int32_t bins, double delta);

  inline int32_t n() const { return n_; }
  inline int32_t p() const { return p_; }
  inline int32_t p2() const { return (p_ - 1) / 2; }
  inline int32_t bins() const { return bins_; }
  inline int32_t *t() const { return t_.data(); }
  inline double *wt() const { return wt_.data(); }

  double SampleInTime(int32_t i) const;
  double SampleInFreq(double xi) const;

private:
  int32_t n_;
  int32_t bins_;
  double delta_;

  double width_;
  double sigma_f_;
  double sigma_t_;
  int32_t p_; // Size of support of window in time domain.

  // We precompute the window in time domain.
  // TODO: It might be more efficient not to store this. Instead, for each t, we
  // can compute the window at t and then iterate over tau's.
  DoubleArray wt_; // Size is p.

  Int32Array t_; // [-p2:p2]. Size is p.
};

} // namespace mps