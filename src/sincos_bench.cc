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

/*
2017-07-22:
Trigonometric operations take up most of the time of binning.
My own sin-cos routine that uses Chebyshev approximation seems fast.
Somehow, the compiler is not SIMD-ing the loops to compute sin-cos.
Hence, my vectorized sin-cos is >10X faster, but it has a ~1e-6 error.
However, Boost-SIMD seems to be the best. It has machine-precision error but
takes the same amount of time.
*/
#include <benchmark/benchmark.h>
#include <eigen3/Eigen/Dense>

#include <boost/align/aligned_allocator.hpp>
#include <boost/simd/arithmetic.hpp>
#include <boost/simd/function/aligned_store.hpp>
#include <boost/simd/function/load.hpp>
#include <boost/simd/function/sin.hpp>
#include <boost/simd/function/sincos.hpp>
#include <boost/simd/function/store.hpp>
#include <boost/simd/memory/allocator.hpp>
#include <boost/simd/pack.hpp>

#include "base.h"

namespace ba = boost::alignment;
namespace bs = boost::simd;

namespace mps {

static void BM_SinTwoPi(benchmark::State &state) {
  const Int n = 4096;
  vector<double> out(n);
  for (int i = 0; i < n; ++i) {
    out[i] = RandomDouble();
  }
  while (state.KeepRunning()) {
    for (int i = 0; i < n; ++i) {
      out[i] = std::sin(out[i] * (2.0 * M_PI));
    }
  }
}
BENCHMARK(BM_SinTwoPi);

static void BM_SinTwoPi_Eigen(benchmark::State &state) {
  const Int n = 4096;
  Eigen::Array<double, n, 1> out;
  for (int i = 0; i < n; ++i) {
    out[i] = RandomDouble();
  }
  Eigen::Array<double, n, 1> out2;
  while (state.KeepRunning()) {
    out2 = out.sin();
  }
}
BENCHMARK(BM_SinTwoPi_Eigen);

static void BM_SinTwoPiApprox(benchmark::State &state) {
  const Int n = 4096;
  DoubleArray out(n);
  for (int i = 0; i < n; ++i) {
    out[i] = RandomDouble();
  }
  double *data = out.data();
  while (state.KeepRunning()) {
    for (int i = 0; i < n; ++i) {
      data[i] = SinTwoPi(data[i]);
    }
  }
}
BENCHMARK(BM_SinTwoPiApprox);

static void BM_SinTwoPiApprox_OMP(benchmark::State &state) {
  const Int n = 4096;
  DoubleArray out(n);
  for (int i = 0; i < n; ++i) {
    out[i] = RandomDouble();
  }
  double *__restrict__ data = out.data();
  while (state.KeepRunning()) {
#pragma omp simd aligned(data : kAlign)
    for (int i = 0; i < n; ++i) {
      data[i] = SinTwoPi(data[i]);
    }
  }
}
BENCHMARK(BM_SinTwoPiApprox_OMP);

static void BM_SinTwoPi_BS(benchmark::State &state) {
  constexpr int n = 4096;
  using pack_t = bs::pack<double>;
  std::vector<double, ba::aligned_allocator<double, pack_t::alignment>> out(n);
  for (int i = 0; i < n; ++i) {
    out[i] = RandomDouble();
  }
  const size_t pack_card = bs::cardinal_of<pack_t>();
  while (state.KeepRunning()) {
    for (size_t i = 0; i < n; i += pack_card) {
      pack_t x(bs::aligned_load<pack_t>(out.data() + i));
      x *= (2 * M_PI);
      x = bs::sin(x);
      bs::aligned_store(x, out.data() + i);
    }
  }
}
BENCHMARK(BM_SinTwoPi_BS);

} // namespace mps