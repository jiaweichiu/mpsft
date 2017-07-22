#include "gen.h"
#include "base.h"
#include "integer.h"
#include "sincos.h"

namespace mps {

void GenerateModeMap(int32_t n, int32_t k, ModeMap *mm) {
  mm->clear();
  const size_t kk = k;
  while (mm->size() < kk) {
    const int32_t idx = PosMod(RandomInt32(), n);
    if (mm->find(idx) != mm->end()) {
      continue;
    }
    Cplex coef(RandomNormal(), RandomNormal());
    coef /= std::abs(coef);
    (*mm)[idx] = coef;
  }
}

void EvaluateModes(int32_t n, const ModeMap &mm, CplexArray *__restrict__ out) {
  CHECK_EQ(out->size(), n);
  out->clear();
  Cplex *data = out->data();

  for (const auto &kv : mm) {
    for (int32_t t = 0; t < n; ++t) {
      const double angle =
          2.0 * M_PI * (double(MulMod(kv.first, t, n)) / double(n));
      data[t] += kv.second * Cplex(::cos(angle), ::sin(angle));
    }
  }
}

// Generate Xhat. Noise energy will sum up to n*sigma*sigma.
void GenerateXhat(int32_t n, const ModeMap &mm, double sigma,
                  CplexArray *__restrict__ out) {
  CHECK_EQ(out->size(), n);
  double noise_energy = 0;
  Cplex *data = out->data();

#pragma omp simd aligned(data : kAlign) reduction(+ : noise_energy)
  for (int i = 0; i < n; ++i) {
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
  for (int i = 0; i < n; ++i) {
    data[i] *= factor;
  }
  for (const auto &kv : mm) {
    (*out)[kv.first] = kv.second;
  }
}

} // namespace mps