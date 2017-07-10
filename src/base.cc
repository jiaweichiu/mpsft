#include <complex>
#include <random>

#include "base.h"

namespace mps {

namespace {

std::mt19937 rng;
std::uniform_int_distribution<Int> uid;

}  // namespace

Int RandomInt() { return uid(rng); }

CplexArray::CplexArray() {}

CplexArray::CplexArray(Int n) { Resize(n); }

void CplexArray::Reset() {
  if (data_) {
    fftw_free(data_);
    n_ = 0;
  }
}

CplexArray::~CplexArray() { Reset(); }

void CplexArray::Resize(Int n) {
  Reset();
  n_ = n;
  data_ = reinterpret_cast<Cplex *>(fftw_malloc(sizeof(Cplex) * n));
}

FFTPlan::FFTPlan(Int n, Int sign, bool in_place) : dummy1_(1), dummy2_(1) {
  fftw_complex *x = reinterpret_cast<fftw_complex *>(dummy1_.data());
  fftw_complex *y = reinterpret_cast<fftw_complex *>(in_place ? dummy1_.data()
                                                          : dummy2_.data());
  plan_ = fftw_plan_dft_1d(n, x, y, sign, FFTW_ESTIMATE);
}

FFTPlan::~FFTPlan() {
  fftw_destroy_plan(plan_);
}

}  // namespace mps