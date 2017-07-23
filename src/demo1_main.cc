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
#include <fstream>
#include <gperftools/profiler.h>

#include "base.h"
#include "demo1.h"

namespace mps {

void Run(int argc, char *const argv[]) {
  MainInit(argc, argv);
  Demo1Options opt;
  opt.trials = 1;
  opt.min_bins = 101;
  opt.window_delta = 1e-5;
  opt.window_threshold = 0.1;
  opt.max_stale_iter = 5;

  constexpr int32_t n = kPrimes[22];
  constexpr int32_t k = 1 << 10;
  constexpr double sigma_min = 1e-3;
  constexpr double sigma_max = 1.0;
  const double l_sigma_min = ::log(sigma_min);
  const double l_sigma_max = ::log(sigma_max);
  constexpr int num_sigma = 6;
  constexpr int num_trials = 10;

  std::ofstream fout("./demo1_main.tsv");
  fout << "i\tj\tsigma\thit\thit_l1_err\thit_l2_err\t"
       << "miss\tmiss_l1_err\tmiss_l2_err\t"
       << "l1_err\tl2_err\tmax_err\n";

  for (int i = 0; i <= num_sigma; ++i) {
    const double l =
        (l_sigma_max - l_sigma_min) * double(i) / double(num_sigma);
    const double sigma = std::exp(l + l_sigma_min);

    for (int j = 0; j < num_trials; ++j) {
      Demo1 demo(opt, n, k, sigma);
      // ProfilerStart("/tmp/demo1_main.prof");
      demo.Run();
      // ProfilerStop();
      fout << i << "\t" << j << "\t" << sigma << "\t";
      demo.PostAnalyze(fout);
    }
  }
}

} // namespace mps

int main(int argc, char *const argv[]) { mps::Run(argc, argv); }