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
#include <complex>
#include <cstdlib>
#include <cstring>
#include <random>

#include "base.h"

namespace mps {

magics_info kPrimesMagic[kNumPrimes];

namespace {

std::mt19937 rng;
std::uniform_int_distribution<int32_t> uid32;
std::uniform_int_distribution<int64_t> uid64;
std::normal_distribution<double> nd;
std::uniform_real_distribution<double> ud;

} // namespace

void MainInit(int argc, char *const argv[]) {
  google::InitGoogleLogging(argv[0]);
  google::InstallFailureSignalHandler();
  rng.seed(123534);
  LOG(INFO) << "mps::MainInit";
  // Compute magic for prime. We will keep dividing by these primes.
  // Don't bother about small primes.
  for (int32_t i = 1; i < kNumPrimes; ++i) {
    kPrimesMagic[i] = compute_signed_magic_info(kPrimes[i]);
  }
}

void RandomSeed(int32_t seed) { rng.seed(seed); }
int32_t RandomInt32() { return uid32(rng); }
int64_t RandomInt64() { return uid64(rng); }
double RandomNormal() { return nd(rng); }
double RandomDouble() { return ud(rng); }

double SincPi(double x) {
  constexpr double taylor_0_bound = std::numeric_limits<double>::epsilon();
  constexpr double taylor_2_bound = std::sqrt(taylor_0_bound);
  constexpr double taylor_n_bound = std::sqrt(taylor_2_bound);
  if (x < 0) {
    x = -x;
  }
  if (x >= taylor_n_bound) {
    return std::sin(x) / x;
  }
  double result = 1.0;
  if (x >= taylor_0_bound) {
    const double x2 = x * x;
    result -= x2 / 6.0;
    if (x >= taylor_2_bound) {
      result += (x2 * x2) / 120.0;
    }
  }
  return result;
}

CplexArray::CplexArray(std::initializer_list<Cplex> l) {
  resize(l.size());
  std::copy(l.begin(), l.end(), data_);
}

double CplexArray::energy() const {
  double ans = 0;
  const Cplex *__restrict__ d = data_;
#pragma omp simd aligned(d : kAlign) reduction(+ : ans)
  for (int32_t i = 0; i < n_; ++i) {
    ans += AbsSq(RE(d[i]), IM(d[i]));
  }
  return ans;
}

CplexMatrix::CplexMatrix(int32_t rows, int32_t cols)
    : rows_(rows), cols_(cols), data_(rows) {
  for (int32_t i = 0; i < rows; ++i) {
    data_[i].resize(cols);
  }
}

void CplexMatrix::clear() {
  for (int32_t i = 0; i < rows_; ++i) {
    data_[i].clear();
  }
}

FFTPlan::FFTPlan(int32_t n, char sign)
    : n_(n), sign_(sign), dummy1_(1), dummy2_(1) {
  fftw_complex *x = reinterpret_cast<fftw_complex *>(dummy1_.data());
  fftw_complex *y = reinterpret_cast<fftw_complex *>(dummy2_.data());
  plan_ = fftw_plan_dft_1d(n, x, y, sign, FFTW_ESTIMATE);
}

FFTPlan::~FFTPlan() { fftw_destroy_plan(plan_); }

void FFTPlan::Run(const CplexArray &u, CplexArray *v) {
  fftw_execute_dft(plan_, reinterpret_cast<fftw_complex *>(u.data()),
                   reinterpret_cast<fftw_complex *>(v->data()));
}

} // namespace mps