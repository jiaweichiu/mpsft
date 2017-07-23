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
#include <iostream>

#include "base.h"

namespace mps {

struct Demo1Options {
  int32_t trials = 1;
  int32_t min_bins = 101;
  double window_delta = 1e-6;
  double window_threshold = 0.1;

  // Exit when we haven't done anything for too many iterations.
  int32_t max_stale_iter = 3;
};

class Demo1 {
public:
  // k is number of modes.
  Demo1(const Demo1Options &opt, int32_t n, int32_t k, double sigma);
  void Run();

  // Reports the errors.
  void PostAnalyze(std::ostream &fout);

private:
  const Demo1Options &opt_;
  int32_t k_;
  CplexArray x_;
  CplexArray xh_;
  double sigma_;

  ModeMap mm_; // Randomly generated.
  ModeMap found_mm_;
};

} // namespace mps