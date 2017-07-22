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
#include <eigen3/Eigen/Dense>

#include "base.h"
#include <benchmark/benchmark.h>

namespace mps {

static void BM_ComplexOp(benchmark::State &state) {
  Cplex x(RandomNormal(), RandomNormal());
  Cplex y(RandomNormal(), RandomNormal());
  Cplex z(0, 0);
  constexpr int n = 2048;
  while (state.KeepRunning()) {
    for (int i = 0; i < n; ++i) {
      z += x * y;
    }
  }
  benchmark::DoNotOptimize(z);
}
BENCHMARK(BM_ComplexOp);

// Evaluate directly. When ffast-math is on, this should have no difference with
// the above.
static void BM_ComplexOp_OMP(benchmark::State &state) {
  Cplex x(RandomNormal(), RandomNormal());
  Cplex y(RandomNormal(), RandomNormal());
  double re = 0;
  double im = 0;
  constexpr int n = 2048;
  while (state.KeepRunning()) {
#pragma omp simd reduction(+ : re, im)
    for (int i = 0; i < n; ++i) {
      re += x.real() * y.real() - x.imag() * y.imag();
      im += x.imag() * y.real() + x.real() * y.imag();
    }
  }
  benchmark::DoNotOptimize(re);
  benchmark::DoNotOptimize(im);
}
BENCHMARK(BM_ComplexOp_OMP);

// Currently, fmod cannot be vectorized.
// No autovectorization expected.
static void BM_PosModOne(benchmark::State &state) {
  const int n = 4096;
  DoubleArray out(n);
  for (int i = 0; i < n; ++i) {
    out[i] = RandomNormal() * 10000;
  }
  double *data = out.data();
  while (state.KeepRunning()) {
    for (int i = 0; i < n; ++i) {
      data[i] = std::fmod(data[i], 1.0);
    }
  }
}
BENCHMARK(BM_PosModOne);

static void BM_PosModOne_OMP(benchmark::State &state) {
  const int n = 4096;
  DoubleArray out(n);
  for (int i = 0; i < n; ++i) {
    out[i] = RandomNormal() * 10000;
  }
  double *__restrict__ data = out.data();
  while (state.KeepRunning()) {
#pragma omp simd aligned(data : kAlign)
    for (int i = 0; i < n; ++i) {
      data[i] = PosModOne(data[i]);
    }
  }
}
BENCHMARK(BM_PosModOne_OMP);

} // namespace mps