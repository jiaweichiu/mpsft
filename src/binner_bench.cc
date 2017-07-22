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
#include "binner.h"
#include "gen.h"
#include "window.h"

namespace mps {

static void BM_BinInTime(benchmark::State &state) {
  const int bin_in_time_type = state.range(0);
  const int32_t n = kPrimes[state.range(1)];
  const int32_t bins = 501;
  const int32_t bits = 15;
  Window win(n, bins, 1e-6);
  CplexMatrix a(1 + 2 * bits, bins);
  std::unique_ptr<BinInTime> bin_in_time(
      BinInTime::Create(bin_in_time_type, win, bits));

  const ModeMap mm = {{5, Cplex(2.0, 1.0)}};
  CplexArray x(n);
  EvaluateModes(n, mm, &x);

  while (state.KeepRunning()) {
    const int32_t q = RandomInt32() % n;
    Transform tf(n);
    bin_in_time->Run(x, tf, q, &a);
  }
}
BENCHMARK(BM_BinInTime)
    ->Args({0, 22})
    ->Args({1, 22})
    ->Args({2, 22})
    ->Args({3, 22})
    ->Args({4, 22});

static void BM_BinInFreq(benchmark::State &state) {
  const int bin_in_freq_type = state.range(0);
  const int32_t n = kPrimes[state.range(1)];
  const int32_t bins = 501;
  const int32_t bits = 15;
  Window win(n, bins, 1e-6);
  CplexMatrix a(1 + 2 * bits, bins);
  std::unique_ptr<BinInFreq> bin_in_freq(
      BinInFreq::Create(bin_in_freq_type, win, bits));

  // Add some modes.
  ModeMap mm;
  GenerateModeMap(n, 100, &mm);

  while (state.KeepRunning()) {
    const int32_t q = RandomInt32() % n;
    Transform tf(n);
    bin_in_freq->Run(mm, tf, q, &a);
  }
}
BENCHMARK(BM_BinInFreq)->Args({0, 22})->Args({1, 22});

} // namespace mps