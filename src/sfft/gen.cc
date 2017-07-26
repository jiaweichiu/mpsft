#include <glog/logging.h>
#include <stdint.h>
#include <stdlib.h>
#include <math.h>

#include "rand.h"
#include "sfft.h"

namespace mps {

void GenerateModeMap(int32_t n, int32_t k, sfft_output *mm) {
  mm->clear();
  const size_t kk = k;
  while (mm->size() < kk) {
    const int32_t idx = ((RandomInt32() % n) + n) % n;
    sfft_output::const_iterator it = mm->find(idx);
    if (it != mm->end()) {
      continue;
    }
    complex_t coef = RandomNormal() + RandomNormal() * I;
    coef /= cabs(coef);
    CHECK(fabs(cabs(coef) - 1) < 1e-4);
    (*mm)[idx] = coef;
  }
}

void GenerateXhatAlt(int32_t n, const sfft_output &mm, double sigma,
                     complex_t *out) {
  memset(out, 0, sizeof(complex_t) * n);
  sfft_output::const_iterator it;
  for (it = mm.begin(); it != mm.end(); ++it) {
    out[it->first] = it->second;
  }
  const double scale = sigma / sqrt(2 * n);
  for (int i = 0; i < n; ++i) {
    complex_t delta = RandomNormal() + RandomNormal() * I;
    out[i] += (delta * scale);
  }
}

}  // namespace mps