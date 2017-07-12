#pragma once

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdint>
#include <fftw3.h>
#include <glog/logging.h>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

namespace mps {

using Int = int32_t;
using Long = int64_t;
using Real = double;
using Cplex = std::complex<Real>;

using IntPair = std::pair<Int, Int>;
using RealPair = std::pair<Real, Real>;
using CplexPair = std::pair<Cplex, Cplex>;

using std::string;
using std::unique_ptr;
using std::vector;

Int RandomInt();

// Our only macros! We try not to.
#define RE std::real
#define IM std::imag

// Force cast into longs.
inline Int PosMod(Long x, Long n) { return ((x % n) + n) % n; }
inline Int Mod(Long x, Long n) { return x % n; }

inline Real Square(Real x) { return x * x; }
inline Int Round(Real x) { return std::round(x); }
inline Cplex Sinusoid(Real theta) {
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

  void Resize(Int n);
  void Reset();
  void Fill(Cplex x);
  void Clear();

  inline Int size() const { return n_; }
  inline Cplex &operator[](Int i) { return data_[i]; }
  inline const Cplex &operator[](Int i) const { return data_[i]; }
  inline Cplex *data() const { return data_; }

private:
  Int n_ = 0;
  Cplex *data_ = nullptr;
};

class CplexMatrix {
public:
  CplexMatrix(Int rows, Int cols);
  inline Int rows() const { return rows_; }
  inline Int cols() const { return cols_; }
  inline CplexArray &operator[](Int i) { return data_[i]; }
  inline const CplexArray &operator[](Int i) const { return data_[i]; }

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