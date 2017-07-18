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
#include "base.h"
#include <benchmark/benchmark.h>

namespace mps {

static void BM_ComplexOp(benchmark::State &state) {
  Cplex x(RandomNormal(), RandomNormal());
  Cplex y(RandomNormal(), RandomNormal());
  Cplex z(0, 0);
  while (state.KeepRunning()) {
    for (int i = 0; i < 2000; ++i) {
      z += x * y;
    }
  }
  benchmark::DoNotOptimize(z);
}
BENCHMARK(BM_ComplexOp);

// Evaluate directly. When ffast-math is on, this should have no difference with
// the above.
static void BM_ComplexOpManual(benchmark::State &state) {
  Cplex x(RandomNormal(), RandomNormal());
  Cplex y(RandomNormal(), RandomNormal());
  double re = 0;
  double im = 0;
  while (state.KeepRunning()) {
    for (int i = 0; i < 2000; ++i) {
      re += x.real() * y.real() - x.imag() * y.imag();
      im += x.imag() * y.real() + x.real() * y.imag();
    }
  }
  benchmark::DoNotOptimize(re);
  benchmark::DoNotOptimize(im);
}
BENCHMARK(BM_ComplexOpManual);

static void BM_Sin(benchmark::State &state) {
  const Int n = 1000;
  while (state.KeepRunning()) {
    for (int i = 0; i <= n; ++i) {
      const double x = double(i) / double(n);
      benchmark::DoNotOptimize(std::sin(x * 2.0 * M_PI));
    }
  }
}
BENCHMARK(BM_Sin);

static void BM_SinTwoPi(benchmark::State &state) {
  const Int n = 1000;
  while (state.KeepRunning()) {
    for (int i = 0; i <= n; ++i) {
      const double x = double(i) / double(n);
      benchmark::DoNotOptimize(SinTwoPi(x));
    }
  }
}
BENCHMARK(BM_SinTwoPi);

static void BM_SinCos(benchmark::State &state) {
  const Int n = 1000;
  double s;
  double c;
  while (state.KeepRunning()) {
    for (int i = 0; i <= n; ++i) {
      const double x = double(i) / double(n) * (2.0 * M_PI);
      ::sincos(x, &s, &c);
    }
  }
}
BENCHMARK(BM_SinCos);

static void BM_SinCosTwoPi(benchmark::State &state) {
  const Int n = 1000;
  while (state.KeepRunning()) {
    for (int i = 0; i <= n; ++i) {
      const double x = double(i) / double(n);
      benchmark::DoNotOptimize(SinCosTwoPi(x));
    }
  }
}
BENCHMARK(BM_SinCosTwoPi);

} // namespace mps