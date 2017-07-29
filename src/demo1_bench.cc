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
#include <benchmark/benchmark.h>

#include "base.h"
#include "demo1.h"
#include "iterate.h"

namespace mps {

// Parameters (n, k): n is size of x and there are k modes.
// TODO: Allow sigma to be varied.
static void BM_Demo1(benchmark::State &state) {
  const int32_t n = state.range(0);
  const int32_t num_modes = state.range(1);
  const double sigma = 1e-1;

  Demo1Options opt;
  opt.trials = 1;
  opt.min_bins = 51;
  opt.window_delta = 1e-5;
  opt.window_threshold = 0.1;
  opt.max_stale_iter = 5;

  Demo1 demo(opt, n, num_modes, sigma);
  while (state.KeepRunning()) {
    demo.Run();
  }
}
BENCHMARK(BM_Demo1)
    ->Args({kPrimes[22], 1 << 6})
    ->Args({kPrimes[22], 1 << 7})
    ->Args({kPrimes[22], 1 << 8})
    ->Args({kPrimes[22], 1 << 9})
    ->Args({kPrimes[22], 1 << 10})
    ->Args({kPrimes[22], 1 << 11})
    ->Args({kPrimes[22], 1 << 12});

BENCHMARK(BM_Demo1)
    ->Args({kPrimes[13], 50})
    ->Args({kPrimes[14], 50})
    ->Args({kPrimes[15], 50})
    ->Args({kPrimes[16], 50})
    ->Args({kPrimes[17], 50})
    ->Args({kPrimes[18], 50})
    ->Args({kPrimes[19], 50})
    ->Args({kPrimes[20], 50})
    ->Args({kPrimes[21], 50})
    ->Args({kPrimes[22], 50})
    ->Args({kPrimes[23], 50})
    ->Args({kPrimes[24], 50})
    ->Args({kPrimes[25], 50})
    ->Args({kPrimes[26], 50});

} // namespace mps