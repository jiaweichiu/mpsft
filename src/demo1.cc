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

#include "demo1.h"
#include "base.h"
#include "gen.h"
#include "iterate.h"

namespace mps {

Demo1::Demo1(const Demo1Options &opt, int32_t n, int32_t k, double sigma)
    : opt_(opt), k_(k), x_(n), xh_(n), sigma_(sigma) {
  GenerateModeMap(n, k, &mm_);
  // GenerateXhat(n, mm_, sigma, &xh_);
  GenerateXhatAlt(n, mm_, sigma, &xh_);
  FFTPlan plan(n, FFTW_BACKWARD);
  plan.Run(xh_, &x_);
}

void Demo1::Run() {
  IterateOptions i_opt;
  i_opt.window_delta = opt_.window_delta;
  i_opt.window_threshold = opt_.window_threshold;
  i_opt.trials = opt_.trials;

  VLOG(2) << "Running Demo1";

  int32_t stale_iter = 0;
  found_mm_.clear();

  for (int32_t i = 0; stale_iter < opt_.max_stale_iter; ++i) {
    i_opt.bins =
        std::max<int32_t>((k_ - found_mm_.size()) * 2 + 1, opt_.min_bins);
    i_opt.bin_threshold = 3.0 * sigma_ / std::sqrt(double(i_opt.bins));
    i_opt.sv_threshold = 1.0 * sigma_ / std::sqrt(double(i_opt.bins));

    const bool res = Iterate(x_, i_opt, &found_mm_);
    VLOG(2) << "iter=" << i << " found=" << found_mm_.size() << " res=" << res;
    if (!res) {
      ++stale_iter;
    }
  }
}

void Demo1::PostAnalyze(std::ostream &fout) {
  const int32_t n = x_.size();

  CplexArray xh2(n);
  xh2.clear();

  int32_t hit = 0;
  int32_t miss = 0;
  double hit_l1_sum = 0;
  double hit_l2_sum = 0;
  double miss_l1_sum = 0;
  double miss_l2_sum = 0;
  for (const auto &kv : found_mm_) {
    xh2[kv.first] = kv.second;
    auto it = mm_.find(kv.first);
    const double err = std::abs(kv.second - xh_[kv.first]);
    if (it == mm_.end()) {
      ++miss;
      miss_l1_sum += err;
      miss_l2_sum += err * err;
    } else {
      ++hit;
      hit_l1_sum += err;
      hit_l2_sum += err * err;
    }
  }
  const double hit_l1_err = hit_l1_sum / hit;
  const double hit_l2_err = std::sqrt(hit_l2_sum / hit);
  fout << hit << "\t" << hit_l1_err << "\t" << hit_l2_err << "\t";
  LOG(INFO) << "hit count=" << hit << " l1_err=" << hit_l1_err
            << " l2_err=" << hit_l2_err;

  const double miss_l1_err = miss_l1_sum / miss;
  const double miss_l2_err = std::sqrt(miss_l2_sum / miss);
  fout << miss << "\t" << miss_l1_err << "\t" << miss_l2_err << "\t";
  LOG(INFO) << "miss count=" << miss << " l1_err=" << miss_l1_err
            << " l2_err=" << miss_l2_err;

  double sum1 = 0;
  double sum2 = 0;
  double max_err = 0;
  for (int32_t i = 0; i < n; ++i) {
    const double err = std::abs(xh2[i] - xh_[i]);
    sum1 += err;
    sum2 += err * err;
    max_err = std::max(max_err, err);
  }
  const double l1_err = sum1 / n;
  const double l2_err = std::sqrt(sum2 / n);
  fout << l1_err << "\t" << l2_err << "\t" << max_err << "\n";
  LOG(INFO) << l1_err << " " << l2_err << " " << max_err;
}

} // namespace mps