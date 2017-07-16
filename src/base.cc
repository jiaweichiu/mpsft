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
#include <cstring>
#include <random>

#include "base.h"

namespace mps {

namespace {

std::mt19937 rng;
std::uniform_int_distribution<Int> uid;
std::normal_distribution<double> nd;

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
}

void RandomSeed(Int seed) { rng.seed(seed); }
Int RandomInt() { return uid(rng); }
double RandomNormal() { return nd(rng); }

CplexArray::CplexArray() {}

CplexArray::CplexArray(Int n) { Resize(n); }

CplexArray::CplexArray(std::initializer_list<Cplex> l) {
  Resize(l.size());
  std::copy(l.begin(), l.end(), data_);
}

void CplexArray::Reset() {
  if (data_) {
    fftw_free(data_);
    n_ = 0;
    data_ = nullptr;
  }
}

double CplexArray::Energy() const {
  double ans = 0;
  for (Int i = 0; i < n_; ++i) {
    ans += AbsSq(data_[i]);
  }
  return ans;
}

CplexArray::~CplexArray() { Reset(); }

void CplexArray::Resize(Int n) {
  Reset();
  n_ = n;
  data_ = reinterpret_cast<Cplex *>(fftw_alloc_complex(n));
}

void CplexArray::Fill(Cplex x) { std::fill(data_, data_ + n_, x); }

void CplexArray::Clear() { std::memset(data_, 0, sizeof(Cplex) * n_); }

CplexArray EvaluateModes(Int n, const ModeMap &mm) {
  CplexArray x(n);
  x.Clear();
  for (const auto &kv : mm) {
    for (Int t = 0; t < n; ++t) {
      const double freq = Mod(kv.first * double(t), n) / double(n);
      x[t] += kv.second * Sinusoid(freq);
    }
  }
  return x;
}

// Generate Xhat. Noise energy will sum up to n*sigma*sigma.
CplexArray GenerateXhat(Int n, const ModeMap &mm, double sigma) {
  CplexArray out(n);
  double noise_energy = 0;
  for (Int i = 0; i < n; ++i) {
    out[i] = Cplex(RandomNormal(), RandomNormal());
    noise_energy += AbsSq(out[i]);
  }
  for (const auto &kv : mm) {
    noise_energy -= AbsSq(out[kv.first]);
  }
  // Rescale noise by this factor.
  const double factor = sigma / std::sqrt(noise_energy);
  for (Int i = 0; i < n; ++i) {
    out[i] *= factor;
  }
  for (const auto &kv : mm) {
    out[kv.first] = kv.second;
  }
  return out;
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

Transform::Transform(Int n) {
  a = (RandomInt() % (n - 1)) + 1;
  b = RandomInt() % n;
  c = RandomInt() % n;
  a_inv = InvMod(a, n);
}

Transform::Transform(Int n, Int a, Int b, Int c) {
  this->a = a;
  this->b = b;
  this->c = c;
  this->a_inv = InvMod(a, n);
}

} // namespace mps