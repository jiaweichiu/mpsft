#include <complex>
#include <random>

#include "base.h"

namespace mps {

namespace {

std::mt19937 rng;
std::uniform_int_distribution<Int> uid;
std::normal_distribution<Real> nd;

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

void MainInit(int argc, char *const argv[]) {
  google::InitGoogleLogging(argv[0]);
  google::InstallFailureSignalHandler();
  rng.seed(123534);
  LOG(INFO) << "mps::MainInit";
}

void RandomSeed(Long seed) { rng.seed(seed); }
Int RandomInt() { return uid(rng); }
Real RandomNormal() { return nd(rng); }

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

CplexArray::~CplexArray() { Reset(); }

void CplexArray::Resize(Int n) {
  Reset();
  n_ = n;
  data_ = reinterpret_cast<Cplex *>(fftw_alloc_complex(n));
}

void CplexArray::Fill(Cplex x) { std::fill(data_, data_ + n_, x); }

void CplexArray::Clear() { Fill(Cplex(0, 0)); }

CplexArray EvaluateModes(Int n, const ModeMap &mm) {
  CplexArray x(n);
  x.Clear();
  for (const auto& kv : mm) {
    for (Int t = 0; t < n; ++t) {
      const Real angle = (2.0 * M_PI) * Mod(kv.first * Real(t), n) / Real(n);
      x[t] += kv.second * Sinusoid(angle);
    }
  }
  return x;
}

// Generate Xhat with desired SNR.
// SNR defined here as sqrt(signal energy / noise energy).
CplexArray GenerateXhat(Int n, const ModeMap& mm, Real snr) {
  CplexArray out(n);
  Real noise_energy = 0;
  for (Int i = 0; i < n; ++i) {
    out[i] = Cplex(RandomNormal(), RandomNormal());
    noise_energy += AbsSq(out[i]);
  }
  Real signal_energy = 0;
  for (const auto& kv : mm) {
    signal_energy += AbsSq(kv.second);
    noise_energy -= AbsSq(out[kv.first]);
  }
  // Rescale noise by this factor.
  const Real factor = std::sqrt(signal_energy / noise_energy) / snr;
  for (Int i = 0; i < n; ++i) {
    out[i] *= factor;
  }
  for (const auto& kv : mm) {
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