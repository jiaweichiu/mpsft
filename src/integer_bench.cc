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

/*
2017-07-22:
AVX2 integer intrinsics only support int32 but not int64.
We are usually dividing/modding by int32 but the product can extend 32-bit.
*/

#include <benchmark/benchmark.h>

#include "base.h"
#include "integer.h"

namespace mps {

const int n = 8192;

static void BM_Int32Div_OMP(benchmark::State &state) {
  Int32Array out(n);
  for (int i = 0; i < n; ++i) {
    out[i] = RandomInt32();
  }
  int32_t *__restrict__ data = out.data();
  while (state.KeepRunning()) {
#pragma omp simd aligned(data : kAlign)
    for (int i = 0; i < n; ++i) {
      data[i] = data[i] / int32_t(kPrimes[15]);
    }
  }
}
BENCHMARK(BM_Int32Div_OMP);

static void BM_Int64Div_OMP(benchmark::State &state) {
  Int64Array out(n);
  for (int i = 0; i < n; ++i) {
    out[i] = RandomInt64();
  }
  int64_t *__restrict__ data = out.data();
  while (state.KeepRunning()) {
#pragma omp simd aligned(data : kAlign)
    for (int i = 0; i < n; ++i) {
      data[i] = data[i] / kPrimes[15];
    }
  }
}
BENCHMARK(BM_Int64Div_OMP);

static void BM_Int32Div_Eigen(benchmark::State &state) {
  Eigen::Array<int32_t, n, 1> out;
  // Int64Array out(n);
  for (int i = 0; i < n; ++i) {
    out[i] = RandomInt32();
  }
  // int64_t * __restrict__ data = out.data();
  while (state.KeepRunning()) {
    out /= kPrimes[15];
  }
}
BENCHMARK(BM_Int32Div_Eigen);

static void BM_Int64Div_Eigen(benchmark::State &state) {
  Eigen::Array<int64_t, n, 1> out;
  // Int64Array out(n);
  for (int i = 0; i < n; ++i) {
    out[i] = RandomInt64();
  }
  // int64_t * __restrict__ data = out.data();
  while (state.KeepRunning()) {
    out /= kPrimes[15];
  }
}
BENCHMARK(BM_Int64Div_Eigen);

static void BM_Mod(benchmark::State &state) {
  const int32_t divisor = kPrimes[15];
  Int32Array out(n);
  for (int i = 0; i < n; ++i) {
    out[i] = RandomInt32();
  }
  int32_t *__restrict__ data = out.data();
  while (state.KeepRunning()) {
#pragma omp simd aligned(data : kAlign)
    for (int i = 0; i < n; ++i) {
      data[i] %= divisor;
    }
  }
}
BENCHMARK(BM_Mod);

static void BM_PosMod1(benchmark::State &state) {
  const int32_t divisor = kPrimes[15];
  Int32Array out(n);
  for (int i = 0; i < n; ++i) {
    out[i] = RandomInt32();
  }
  int32_t *__restrict__ data = out.data();
  while (state.KeepRunning()) {
#pragma omp simd aligned(data : kAlign)
    for (int i = 0; i < n; ++i) {
      data[i] = UnsafePosMod1(data[i], divisor);
    }
  }
}
BENCHMARK(BM_PosMod1);

static void BM_PosMod2(benchmark::State &state) {
  const int32_t divisor = kPrimes[15];
  Int32Array out(n);
  for (int i = 0; i < n; ++i) {
    out[i] = RandomInt32();
  }
  int32_t *__restrict__ data = out.data();
  while (state.KeepRunning()) {
#pragma omp simd aligned(data : kAlign)
    for (int i = 0; i < n; ++i) {
      data[i] = UnsafePosMod2(data[i], divisor);
    }
  }
}
BENCHMARK(BM_PosMod2);

static void BM_PosMod(benchmark::State &state) {
  const int32_t divisor = kPrimes[15];
  Int32Array out(n);
  for (int i = 0; i < n; ++i) {
    out[i] = RandomInt32();
  }
  int32_t *__restrict__ data = out.data();
  while (state.KeepRunning()) {
#pragma omp simd aligned(data : kAlign)
    for (int i = 0; i < n; ++i) {
      data[i] = UnsafePosMod(data[i], divisor);
    }
  }
}
BENCHMARK(BM_PosMod);

} // namespace mps