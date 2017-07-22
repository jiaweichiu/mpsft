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
#pragma once

// "Base" lib contains a bunch of handy math functions.

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdint>
#include <fftw3.h>
#include <glog/logging.h>
#include <initializer_list>
#include <iostream>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

#include "magic.h"

namespace mps {

constexpr size_t kAlign = 64;
using Cplex = std::complex<double>;

using ModeMap = std::unordered_map<int32_t, Cplex>;

using std::string;
using std::vector;

// List of primes close to powers of 2.
constexpr int kNumPrimes = 30;
constexpr int32_t kPrimes[kNumPrimes] = {
    2,        3,        5,        7,         17,        31,
    67,       127,      257,      509,       1021,      2053,
    4099,     8191,     16381,    32771,     65537,     131071,
    262147,   524287,   1048573,  2097143,   4194301,   8388617,
    16777213, 33554467, 67108859, 134217757, 268435459, 536870909};

extern magics_info kPrimesMagic[kNumPrimes];

void MainInit(int argc, char *const argv[]);

void RandomSeed(int32_t seed);
int32_t RandomInt32();
int64_t RandomInt64();
double RandomDouble(); // [0, 1).
double RandomNormal();

// Our only macros! We try not to.
#define RE std::real
#define IM std::imag

#pragma omp declare simd
inline double PosModOne(double x) { return x - std::floor(x); }

// #pragma omp declare simd
// inline double AbsSq(double re, double im) { return re * re + im * im; }

// #pragma omp declare simd
// inline double Square(double x) { return x * x; }

// Equivalent to multiplying by i: (x+iy)*i = -y+ix.
// inline Cplex RotateForward(Cplex x) { return Cplex(-IM(x), RE(x)); }

// Equivalent to multiplying by -i: (x+iy)*(-i) = y-ix.
// inline Cplex RotateBackward(Cplex x) { return Cplex(IM(x), -RE(x)); }

template <class T> class Array {
public:
  Array(int n) : n_(n) {
    size_t m = n * sizeof(T);
    m = ((m + kAlign - 1) / kAlign) * kAlign;
    data_ = reinterpret_cast<T *>(::aligned_alloc(kAlign, m));
    CHECK(data_);
  }
  ~Array() { ::free(data_); }

  inline int size() const { return n_; }
  inline T &operator[](int i) { return data_[i]; }
  inline const double &operator[](int i) const { return data_[i]; }
  inline T *data() const { return data_; }

private:
  int n_ = 0;
  T *data_ = nullptr;
};
using DoubleArray = Array<double>;
using Int32Array = Array<int32_t>;
using Int64Array = Array<int64_t>;

class CplexArray {
public:
  CplexArray() = default;
  CplexArray(int32_t n);
  ~CplexArray();
  CplexArray(std::initializer_list<Cplex> l);

  void resize(int32_t n);
  void reset();
  void fill(Cplex x);
  void clear();
  double energy() const;

  inline int32_t size() const { return n_; }
  inline Cplex &operator[](int32_t i) { return data_[i]; }
  inline const Cplex &operator[](int32_t i) const { return data_[i]; }
  inline Cplex *data() const { return data_; }

private:
  int32_t n_ = 0;
  Cplex *data_ = nullptr;
};

class CplexMatrix {
public:
  CplexMatrix(int32_t rows, int32_t cols);
  inline int32_t rows() const { return rows_; }
  inline int32_t cols() const { return cols_; }
  inline CplexArray &operator[](int32_t i) { return data_[i]; }
  inline const CplexArray &operator[](int32_t i) const { return data_[i]; }
  void clear();

private:
  int32_t rows_;
  int32_t cols_;
  vector<CplexArray> data_;
};

class FFTPlan {
public:
  // For sign:
  // #define FFTW_FORWARD (-1)
  // #define FFTW_BACKWARD (+1)
  FFTPlan(int32_t n, char sign);
  ~FFTPlan();

  inline int32_t n() const { return n_; }
  inline char sign() const { return sign_; }

  // If in-place, u, v should be the same.
  // We assume this agrees with in_place supplied to constructor.
  void Run(const CplexArray &u, CplexArray *v);

private:
  int32_t n_;
  char sign_;
  CplexArray dummy1_, dummy2_;
  fftw_plan plan_;
};

} // namespace mps