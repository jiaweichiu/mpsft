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
  Window(Int n, Int bins, double delta);

  // Assume 0 <= t <= p2 where p2=(p-1)/2. Assume p is odd.
  inline double wt(Int t) const {
    DCHECK_GE(t, 0);
    DCHECK_LE(t, p2());
    return wt_[t];
  }

  inline Int n() const { return n_; }
  inline Int p() const { return p_; }
  inline Int p2() const { return (p_ - 1) / 2; }
  inline Int bins() const { return bins_; }

  // double Energy() const;

  double SampleInTime(Int i) const;
  double SampleInFreq(double xi) const;

private:
  Int n_;
  Int bins_;
  double delta_;

  double width_;
  double sigma_f_;
  double sigma_t_;
  Int p_; // Size of support of window in time domain.

  // We precompute the window in time domain.
  // TODO: It might be more efficient not to store this. Instead, for each t, we
  // can compute the window at t and then iterate over tau's.
  vector<double> wt_; // Size is (p-1)/2.
};

} // namespace mps