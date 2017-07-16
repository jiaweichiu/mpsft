#pragma once

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

namespace mps {

// using Int = int32_t;
using Int = int64_t;
using Long = int64_t;
using Cplex = std::complex<double>;

using ModeMap = std::unordered_map<Int, Cplex>;

using std::string;
using std::unique_ptr;
using std::vector;

// List of primes close to powers of 2.
constexpr Int kPrimes[] = {2,        3,        5,         7,         17,
                           31,       67,       127,       257,       509,
                           1021,     2053,     4099,      8191,      16381,
                           32771,    65537,    131071,    262147,    524287,
                           1048573,  2097143,  4194301,   8388617,   16777213,
                           33554467, 67108859, 134217757, 268435459, 536870909};

void MainInit(int argc, char *const argv[]);

void RandomSeed(Long seed);
Int RandomInt();
double RandomNormal();

// Our only macros! We try not to.
#define RE std::real
#define IM std::imag

// Force cast into longs.
inline Int PosMod(Long x, Long n) { return ((x % n) + n) % n; }
inline Int Mod(Long x, Long n) { return x % n; }
inline double AbsSq(Cplex x) { return RE(x) * RE(x) + IM(x) * IM(x); }

inline double Square(double x) { return x * x; }
inline Cplex Sinusoid(double theta) {
  return Cplex(std::cos(theta), std::sin(theta));
}

// Equivalent to multiplying by i: (x+iy)*i = -y+ix.
inline Cplex RotateForward(Cplex x) { return Cplex(-IM(x), RE(x)); }

// Equivalent to multiplying by -i: (x+iy)*(-i) = y-ix.
inline Cplex RotateBackward(Cplex x) { return Cplex(IM(x), -RE(x)); }

class CplexArray {
public:
  CplexArray();
  CplexArray(Int n);
  ~CplexArray();
  CplexArray(std::initializer_list<Cplex> l);

  void Resize(Int n);
  void Reset();
  void Fill(Cplex x);
  void Clear();
  double Energy() const;

  inline Int size() const { return n_; }
  inline Cplex &operator[](Int i) { return data_[i]; }
  inline const Cplex &operator[](Int i) const { return data_[i]; }
  inline Cplex *data() const { return data_; }

private:
  Int n_ = 0;
  Cplex *data_ = nullptr;
};

CplexArray EvaluateModes(Int n, const ModeMap &mm);

// Add ambience noise such that in the *time domain*, each sample point is
// contaminated by N(0, sigma).
// Note: x(t) = sum_k xh[k] exp(2*pi*i*k*t). This is unnormalized.
// If xh[k] ~ N(0, s*s), then x(t) ~ N(0, s*s*n) where s*s*n=sigma*sigma.
CplexArray GenerateXhat(Int n, const ModeMap &mm, double sigma);

class CplexMatrix {
public:
  CplexMatrix(Int rows, Int cols);
  inline Int rows() const { return rows_; }
  inline Int cols() const { return cols_; }
  inline CplexArray &operator[](Int i) { return data_[i]; }
  inline const CplexArray &operator[](Int i) const { return data_[i]; }
  void Clear();

private:
  Int rows_;
  Int cols_;
  vector<CplexArray> data_;
};

class FFTPlan {
public:
  // For sign:
  // #define FFTW_FORWARD (-1)
  // #define FFTW_BACKWARD (+1)
  FFTPlan(Int n, char sign);
  ~FFTPlan();

  inline Int n() const { return n_; }
  inline char sign() const { return sign_; }

  // If in-place, u, v should be the same.
  // We assume this agrees with in_place supplied to constructor.
  void Run(const CplexArray &u, CplexArray *v);

private:
  Int n_;
  char sign_;
  CplexArray dummy1_, dummy2_;
  fftw_plan plan_;
};

// y[t] = x[a*t+c] exp(2*pi*i*b*t/n).
// yh[a*k+b] = xh[k] exp(2*pi*i*c*k/n).
// Mode permutation: a, b
// Mode modulation: c
// Forward: a*k+b.
// Backward: a_inv*(k-b).
// Very lightweight class. Keep logic in the other code.
struct Transform {
  Transform(Int n);
  Transform(Int n, Int a, Int b, Int c);
  Int a;
  Int b;
  Int c;
  Int a_inv;
};

} // namespace mps