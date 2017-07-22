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
2017-07-21:
Trigonometric operations take up most of the time of binning.
My own sin-cos routine that uses Chebyshev approximation seems fast.
Somehow, the compiler is not SIMD-ing the loops to compute sin-cos.
Hence, my vectorized sin-cos is >10X faster, but it has a ~1e-6 error.
However, Boost-SIMD seems to be the best. It has machine-precision error but
takes the same amount of time.

2017-07-22:
Prefer working with CplexArrays. With that constraint, Boost-SIMD requires
splicing of results back and becomes slightly slower. Somehow, gcc can vectorize
the approx functions properly, without having to splice. Code will look simpler.
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
#include "sincos.h"

namespace ba = boost::alignment;
namespace bs = boost::simd;

namespace mps {

constexpr int n = 8192;

static void BM_SinTwoPi(benchmark::State &state) {
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
  DoubleArray out(n);
  for (int i = 0; i < n; ++i) {
    out[i] = RandomDouble();
  }
  double *data = out.data();
  while (state.KeepRunning()) {
    for (int i = 0; i < n; ++i) {
      data[i] = SinTwoPiApprox(data[i]);
    }
  }
}
BENCHMARK(BM_SinTwoPiApprox);

// Add some OMP pragmas to encourage vectorization.
static void BM_SinTwoPiApprox_OMP(benchmark::State &state) {
  DoubleArray out(n);
  for (int i = 0; i < n; ++i) {
    out[i] = RandomDouble();
  }
  double *__restrict__ data = out.data();
  while (state.KeepRunning()) {
#pragma omp simd aligned(data : kAlign)
    for (int i = 0; i < n; ++i) {
      data[i] = SinTwoPiApprox(data[i]);
    }
  }
}
BENCHMARK(BM_SinTwoPiApprox_OMP);

static void BM_SinTwoPi_BS(benchmark::State &state) {
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

static void BM_SinCosTwoPi(benchmark::State &state) {
  DoubleArray in(n);
  CplexArray out(n);
  for (int i = 0; i < n; ++i) {
    in[i] = RandomDouble();
  }
  while (state.KeepRunning()) {
    for (int i = 0; i < n; ++i) {
      out[i] =
          Cplex(std::cos(in[i] * (2.0 * M_PI)), std::sin(in[i] * (2.0 * M_PI)));
    }
  }
}
BENCHMARK(BM_SinCosTwoPi);

// OMP SIMD. Two separate output arrays.
static void BM_SinCosTwoPiApprox_OMP(benchmark::State &state) {
  DoubleArray in(n);
  DoubleArray a_out1(n);
  DoubleArray a_out2(n);
  for (int i = 0; i < n; ++i) {
    in[i] = RandomDouble();
  }
  double *__restrict__ data = in.data();
  double *__restrict__ out1 = a_out1.data();
  double *__restrict__ out2 = a_out2.data();
  while (state.KeepRunning()) {
#pragma omp simd aligned(data : kAlign)
    for (int i = 0; i < n; ++i) {
      out1[i] = SinTwoPiApprox(data[i]);
      out2[i] = CosTwoPiApprox(data[i]);
    }
  }
}
BENCHMARK(BM_SinCosTwoPiApprox_OMP);

// OMP SIMD. One single output array.
static void BM_SinCosTwoPiApproxCplex_OMP(benchmark::State &state) {
  DoubleArray in(n);
  CplexArray a_out(n);
  for (int i = 0; i < n; ++i) {
    in[i] = RandomDouble();
  }
  double *__restrict__ data = in.data();
  Cplex *__restrict__ out = a_out.data();
  while (state.KeepRunning()) {
#pragma omp simd aligned(data : kAlign)
    for (int i = 0; i < n; ++i) {
      out[i] = Cplex(CosTwoPiApprox(data[i]), SinTwoPiApprox(data[i]));
    }
  }
}
BENCHMARK(BM_SinCosTwoPiApproxCplex_OMP);

// Boost SIMD sincos is very fast but it produces two arrays instead of one
// CplexArray. Let's also measure time to merge them back.
static void BM_SinCosTwoPi_BS(benchmark::State &state) {
  const int splice_back = state.range(0);
  using pack_t = bs::pack<double>;
  std::vector<double, ba::aligned_allocator<double, pack_t::alignment>> in(n);
  std::vector<double, ba::aligned_allocator<double, pack_t::alignment>> out1(n);
  std::vector<double, ba::aligned_allocator<double, pack_t::alignment>> out2(n);
  CplexArray out(n);
  for (int i = 0; i < n; ++i) {
    in[i] = RandomDouble();
  }
  const size_t pack_card = bs::cardinal_of<pack_t>();
  while (state.KeepRunning()) {
    for (size_t i = 0; i < n; i += pack_card) {
      pack_t x(bs::aligned_load<pack_t>(in.data() + i));
      auto res = bs::sincos(2.0 * M_PI * x);
      bs::aligned_store(res.first, out1.data() + i);
      bs::aligned_store(res.second, out2.data() + i);
    }
    if (splice_back) {
      for (size_t i = 0; i < n; ++i) {
        out[i] = Cplex(out1[i], out2[i]);
      }
    }
  }
}
BENCHMARK(BM_SinCosTwoPi_BS)->Arg(0)->Arg(1);

} // namespace mps