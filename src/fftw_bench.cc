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

namespace mps {

static void BM_FFTW(benchmark::State &state) {
  const Int n = state.range(0);
  const Int num_modes = 100; // Shouldn't matter.
  const double sigma = 0.1;

  ModeMap mm;
  GenerateModeMap(n, num_modes, &mm);

  CplexArray xh(n);
  GenerateXhat(n, mm, sigma, &xh);
  CplexArray x(n);

  FFTPlan plan(n, FFTW_BACKWARD);
  while (state.KeepRunning()) {
    plan.Run(xh, &x);
  }
}
BENCHMARK(BM_FFTW)->RangeMultiplier(2)->Range(1 << 9, 1 << 22);

BENCHMARK(BM_FFTW)
    ->Arg(kPrimes[9])
    ->Arg(kPrimes[10])
    ->Arg(kPrimes[11])
    ->Arg(kPrimes[12])
    ->Arg(kPrimes[13])
    ->Arg(kPrimes[14])
    ->Arg(kPrimes[15])
    ->Arg(kPrimes[16])
    ->Arg(kPrimes[17])
    ->Arg(kPrimes[18])
    ->Arg(kPrimes[19])
    ->Arg(kPrimes[20])
    ->Arg(kPrimes[21])
    ->Arg(kPrimes[22]);

} // namespace mps