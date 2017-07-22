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

// PowMod returns mod(b^e, m).
// b=base, e=exponent, m=modulus.
Int PowMod(Int b, Int e, Int m) {
  Int r = 1;
  while (e > 0) {
    if (e & 1) { // Odd exponent
      r = Mod(r * b, m);
    }
    e >>= 1;
    b = Mod(b * b, m);
  }
  return r;
}

// Returns a's inverse modulo m. Caution: Assumes prime m.
Int InvMod(Int a, Int m) { return PowMod(a, m - 2, m); }

} // namespace

void MainInit(int argc, char *const argv[]) {
  google::InitGoogleLogging(argv[0]);
  google::InstallFailureSignalHandler();
  rng.seed(123534);
  LOG(INFO) << "mps::MainInit";
  // Compute magic for prime. We will keep dividing by these primes.
  // Don't bother about small primes.
  for (Int i = 1; i < kNumPrimes; ++i) {
    kPrimesMagic[i] = compute_signed_magic_info(kPrimes[i]);
  }
}

void RandomSeed(Int seed) { rng.seed(seed); }
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

double SinTwoPi(double a) {
  const double y = (a >= 0.75) ? -1.0 + a : (a >= 0.25) ? 0.5 - a : a;
  const double x = 4.0 * y;
  constexpr double coef1 = 1.1336481778117871;
  constexpr double coef3 = -0.13807177660048911;
  constexpr double coef5 = 0.0044907175846143066;
  constexpr double coef7 = -6.8290405376023045e-05;
  constexpr double z0 = 1.0;
  const double z1 = x;
  const double z2 = 2.0 * x * z1 - z0;
  const double z3 = 2.0 * x * z2 - z1;
  const double z4 = 2.0 * x * z3 - z2;
  const double z5 = 2.0 * x * z4 - z3;
  const double z6 = 2.0 * x * z5 - z4;
  const double z7 = 2.0 * x * z6 - z5;
  return coef1 * z1 + coef3 * z3 + coef5 * z5 + coef7 * z7;
}

double CosTwoPi(double x) {
  if (x >= 0.75) {
    x = -1.0 + x;
  } else if (x >= 0.25) {
    x = 0.5 - x;
  }
  x *= 4.0;
  constexpr double coef1 = 1.1336481778117871;
  constexpr double coef3 = -0.13807177660048911;
  constexpr double coef5 = 0.0044907175846143066;
  constexpr double coef7 = -6.8290405376023045e-05;
  constexpr double z0 = 1.0;
  const double z1 = x;
  const double z2 = 2.0 * x * z1 - z0;
  const double z3 = 2.0 * x * z2 - z1;
  const double z4 = 2.0 * x * z3 - z2;
  const double z5 = 2.0 * x * z4 - z3;
  const double z6 = 2.0 * x * z5 - z4;
  const double z7 = 2.0 * x * z6 - z5;
  return coef1 * z1 + coef3 * z3 + coef5 * z5 + coef7 * z7;
}

// Use cos(x)=sin(pi/2+x). Basically, shift by one quadrant.
/*Cplex Sinusoid(double x) {
  if (x >= 1.0) {
    x -= 1.0;
  } else if (x < -1.0) {
    x += 2.0;
  } else if (x < 0) {
    x += 1.0;
  }
  DCHECK_GE(x, 0);
  DCHECK_LE(x, 1.0);

  double xb;
  if (x >= 0.75) {
    // Fourth quadrant.
    xb = x - 0.75; // x+0.25-1.0.
    x = -1.0 + x;  // -0.25 < x < 0.
  } else if (x >= 0.5) {
    // Third quadrant.
    xb = x - 0.75; // -1.0+(x+0.25)=x-0.75.
    x = 0.5 - x;   // -0.25 < x < 0.
  } else if (x >= 0.25) {
    // Second quadrant.
    xb = 0.25 - x; // 0.5-(x+0.25)=0.25-x.
    x = 0.5 - x;   // 0 < x < 0.25.
  } else {
    // First quadrant.
    // 0 < x < 0.25.
    xb = 0.25 - x; // 0.5-(x+0.25)=0.25-x.
  }
  x *= 4.0;
  xb *= 4.0;
  constexpr double coef1 = 1.1336481778117871;
  constexpr double coef3 = -0.13807177660048911;
  constexpr double coef5 = 0.0044907175846143066;
  constexpr double coef7 = -6.8290405376023045e-05;

  const double z1 = x;
  const double z1b = xb;

  const double z2 = 2.0 * x * z1 - 1.0;
  const double z2b = 2.0 * xb * z1b - 1.0;

  const double z3 = 2.0 * x * z2 - z1;
  const double z3b = 2.0 * xb * z2b - z1b;

  const double z4 = 2.0 * x * z3 - z2;
  const double z4b = 2.0 * xb * z3b - z2b;

  const double z5 = 2.0 * x * z4 - z3;
  const double z5b = 2.0 * xb * z4b - z3b;

  const double z6 = 2.0 * x * z5 - z4;
  const double z6b = 2.0 * xb * z5b - z4b;

  const double z7 = 2.0 * x * z6 - z5;
  const double z7b = 2.0 * xb * z6b - z5b;

  const double s = coef1 * z1 + coef3 * z3 + coef5 * z5 + coef7 * z7;
  const double c = coef1 * z1b + coef3 * z3b + coef5 * z5b + coef7 * z7b;
  return Cplex(c, s);
}*/

CplexArray::CplexArray() {}

CplexArray::CplexArray(Int n) { Resize(n); }

CplexArray::CplexArray(std::initializer_list<Cplex> l) {
  Resize(l.size());
  std::copy(l.begin(), l.end(), data_);
}

void CplexArray::Reset() {
  if (data_) {
    ::free(data_);
    n_ = 0;
    data_ = nullptr;
  }
}

double CplexArray::Energy() const {
  double ans = 0;
  Cplex *__restrict__ x = data_;
#pragma omp simd aligned(x : kAlign) reduction(+ : ans)
  for (Int i = 0; i < n_; ++i) {
    ans += AbsSq(RE(x[i]), IM(x[i]));
  }
  return ans;
}

CplexArray::~CplexArray() { Reset(); }

void CplexArray::Resize(Int n) {
  Reset();
  n_ = n;
  size_t m = n * sizeof(Cplex);
  m = ((m + kAlign - 1) / kAlign) * kAlign;
  data_ = reinterpret_cast<Cplex *>(::aligned_alloc(kAlign, m));
  CHECK(data_);
}

void CplexArray::Fill(Cplex x) { std::fill(data_, data_ + n_, x); }

void CplexArray::Clear() { std::memset(data_, 0, sizeof(Cplex) * n_); }

void GenerateModeMap(Int n, Int k, ModeMap *mm) {
  mm->clear();
  const size_t kk = k;
  while (mm->size() < kk) {
    const Int idx = PosMod(RandomInt32(), n);
    if (mm->find(idx) != mm->end()) {
      continue;
    }
    Cplex coef(RandomNormal(), RandomNormal());
    coef /= std::abs(coef);
    (*mm)[idx] = coef;
  }
}

void EvaluateModes(Int n, const ModeMap &mm, CplexArray *__restrict__ out) {
  CHECK_EQ(out->size(), n);
  out->Clear();
  Cplex *data = out->data();

  for (const auto &kv : mm) {
#pragma omp simd aligned(data : kAlign)
    for (Int t = 0; t < n; ++t) {
      const double freq = Mod(kv.first * double(t), n) / double(n);
      data[t] += kv.second * Sinusoid(freq);
    }
  }
}

// Generate Xhat. Noise energy will sum up to n*sigma*sigma.
void GenerateXhat(Int n, const ModeMap &mm, double sigma,
                  CplexArray *__restrict__ out) {
  CHECK_EQ(out->size(), n);
  double noise_energy = 0;
  Cplex *data = out->data();

#pragma omp simd aligned(data : kAlign) reduction(+ : noise_energy)
  for (Int i = 0; i < n; ++i) {
    data[i] = Cplex(RandomNormal(), RandomNormal());
    noise_energy += AbsSq(RE(data[i]), IM(data[i]));
  }
  for (const auto &kv : mm) {
    Cplex x = (*out)[kv.first];
    noise_energy -= AbsSq(RE(x), IM(x));
  }
  // Rescale noise by this factor.
  const double factor = sigma / std::sqrt(noise_energy);

#pragma omp simd aligned(data : kAlign)
  for (Int i = 0; i < n; ++i) {
    data[i] *= factor;
  }
  for (const auto &kv : mm) {
    (*out)[kv.first] = kv.second;
  }
}

CplexMatrix::CplexMatrix(Int rows, Int cols)
    : rows_(rows), cols_(cols), data_(rows) {
  for (Int i = 0; i < rows; ++i) {
    data_[i].Resize(cols);
  }
}

void CplexMatrix::Clear() {
  for (Int i = 0; i < rows_; ++i) {
    data_[i].Clear();
  }
}

FFTPlan::FFTPlan(Int n, char sign)
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

std::ostream &operator<<(std::ostream &dest, __int128_t value) {
  std::ostream::sentry s(dest);
  if (s) {
    __uint128_t tmp = value < 0 ? -value : value;
    char buffer[128];
    char *d = std::end(buffer);
    do {
      --d;
      *d = "0123456789"[tmp % 10];
      tmp /= 10;
    } while (tmp != 0);
    if (value < 0) {
      --d;
      *d = '-';
    }
    int len = std::end(buffer) - d;
    if (dest.rdbuf()->sputn(d, len) != len) {
      dest.setstate(std::ios_base::badbit);
    }
  }
  return dest;
}

} // namespace mps