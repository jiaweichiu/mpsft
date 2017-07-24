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
#include "gen.h"

namespace mps {

static void BM_FFTW(benchmark::State &state) {
  const int32_t n = state.range(0);
  const int flags = state.range(1);
  const int32_t num_modes = 100; // Shouldn't matter.
  const double sigma = 0.1;

  ModeMap mm;
  GenerateModeMap(n, num_modes, &mm);

  CplexArray xh(n);
  GenerateXhat(n, mm, sigma, &xh);
  CplexArray x(n);

  FFTPlan plan(n, FFTW_BACKWARD, flags);
  while (state.KeepRunning()) {
    plan.Run(xh, &x);
  }
}

static void CustomArguments(benchmark::internal::Benchmark *b) {
  std::vector<int> flags = {FFTW_ESTIMATE, FFTW_MEASURE};
  for (int flag : flags) {
    for (int i = 9; i <= 26; ++i) {
      b->Args({1 << i, flag});
    }
  }
}
BENCHMARK(BM_FFTW)->Apply(CustomArguments);

// Don't worry about those primes. We know FFTW is a lot slower on primes.

} // namespace mps