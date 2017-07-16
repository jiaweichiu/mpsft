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
#include "window.h"

namespace mps {

static void BM_BinInTime(benchmark::State &state) {
  const int binner_type = state.range(0);
  const Int n = kPrimes[state.range(1)];
  const Int bins = 5;
  const Int bits = 15;
  Window win(n, bins, 1e-6);
  CplexMatrix a(1 + 2 * bits, bins);
  unique_ptr<Binner> binner(Binner::Create(binner_type, win, bits));

  const ModeMap mm = {{5, Cplex(2.0, 1.0)}};
  const CplexArray x = EvaluateModes(n, mm);

  while (state.KeepRunning()) {
    const Int q = RandomInt() % n;
    Transform tf(n);
    binner->BinInTime(x, tf, q, &a);
  }
}
BENCHMARK(BM_BinInTime)
    ->Args({kBinnerSimple, 10})
    ->Args({kBinnerSimple, 13})
    ->Args({kBinnerSimple, 16})
    ->Args({kBinnerFast, 10})
    ->Args({kBinnerFast, 13})
    ->Args({kBinnerFast, 16});

static void BM_BinInFreq(benchmark::State &state) {
  const int binner_type = state.range(0);
  const Int n = kPrimes[state.range(1)];
  const Int bins = 5;
  const Int bits = 15;
  Window win(n, bins, 1e-6);
  CplexMatrix a(1 + 2 * bits, bins);
  unique_ptr<Binner> binner(Binner::Create(binner_type, win, bits));

  ModeMap mm;
  // Add some modes.
  for (Int i = 0; i < 100; ++i) {
    mm[PosMod(RandomInt(), n)] = Cplex(RandomNormal(), RandomNormal());
  }

  while (state.KeepRunning()) {
    const Int q = RandomInt() % n;
    Transform tf(n);
    binner->BinInFreq(mm, tf, q, &a);
  }
}
BENCHMARK(BM_BinInFreq)
    ->Args({kBinnerSimple, 10})
    ->Args({kBinnerSimple, 13})
    ->Args({kBinnerSimple, 16})
    ->Args({kBinnerFast, 10})
    ->Args({kBinnerFast, 13})
    ->Args({kBinnerFast, 16});

} // namespace mps