#include <complex>
#include <random>

#include "base.h"

namespace mps {

namespace {

std::mt19937 rng;
std::uniform_int_distribution<Int> uid;

// PowMod returns mod(b^e, m).
// b=base, e=exponent, m=modulus.
Int PowMod(Int b, Int e, Int m) {
  Int r = 1;
  while (e > 0) {
    if (e & 1) { // Odd exponent
      r = Mod(Long(r) * Long(b), m);
    }
    e >>= 1;
    b = Mod(Long(b) * Long(b), m);
  }
  return r;
}

// Returns a's inverse modulo m. Caution: Assumes prime m.
Int InvMod(Int a, Int m) { return PowMod(a, m - 2, m); }

} // namespace

Int RandomInt() { return uid(rng); }

CplexArray::CplexArray() {}

CplexArray::CplexArray(Int n) { Resize(n); }

void CplexArray::Reset() {
  if (data_) {
    fftw_free(data_);
    n_ = 0;
    data_ = nullptr;
  }
}

CplexArray::~CplexArray() { Reset(); }

void CplexArray::Resize(Int n) {
  Reset();
  n_ = n;
  data_ = reinterpret_cast<Cplex *>(fftw_alloc_complex(n));
}

void CplexArray::Fill(Cplex x) { std::fill(data_, data_ + n_, x); }

void CplexArray::Clear() { Fill(Cplex(0, 0)); }

FFTPlan::FFTPlan(Int n, Int sign) : n_(n), dummy1_(1), dummy2_(1) {
  fftw_complex *x = reinterpret_cast<fftw_complex *>(dummy1_.Data());
  fftw_complex *y = reinterpret_cast<fftw_complex *>(dummy2_.Data());
  plan_ = fftw_plan_dft_1d(n, x, y, sign, FFTW_ESTIMATE);
}

FFTPlan::~FFTPlan() { fftw_destroy_plan(plan_); }

void FFTPlan::Run(const CplexArray &u, CplexArray *v) {
  fftw_execute_dft(plan_, reinterpret_cast<fftw_complex *>(u.Data()),
                   reinterpret_cast<fftw_complex *>(v->Data()));
}

Transform::Transform(Int n) {
  a = (RandomInt() % (n - 1)) + 1;
  b = RandomInt() % n;
  a_inv = InvMod(a, n);
}

Transform::Transform(Int n, Int a, Int b, Int c) {
  this->a = a;
  this->b = b;
  this->c = c;
  this->a_inv = InvMod(a, n);
}

} // namespace mps