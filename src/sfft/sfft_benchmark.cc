/*
Not using std=c++11 as it does not accept "double complex = complex_t" used in
SFFT3 lib. This is painful as many C++11 features will not available, such as
initializer lists.

tr1 rng classes seem to be buggy. We use Boost instead.
*/
#include <benchmark/benchmark.h>
#include <stdint.h>
#include <complex.h>
#include <math.h>
#include <fftw3.h>
#include <glog/logging.h>
// #include <tr1/random>
#include <vector>

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
    (*mm)[idx] = coef;
  }
}

void GenerateXhatAlt(int32_t n, const sfft_output &mm, double sigma,
                     complex_t *out) {
  memset(out, 0, sizeof(complex_t) * n);
  sfft_output::const_iterator it = mm.begin();
  for (it = mm.begin(); it != mm.end(); ++it) {
    out[it->first] = it->second;
  }
  const double scale = sigma / sqrt(2 * n);
  for (int i = 0; i < n; ++i) {
    complex_t delta = RandomNormal() + RandomNormal() * I;
    out[i] += (delta * scale);
  }
}

static void BM_SFFT(benchmark::State &state) {
  const int32_t n = state.range(0);
  const int32_t k = state.range(1); // Sparsity.
  const sfft_version version = static_cast<sfft_version>(state.range(2));
  const double sigma = 1e-1;

  complex_t *x = (complex_t *)sfft_malloc(sizeof(complex_t) * n);
  complex_t *xh = (complex_t *)sfft_malloc(sizeof(complex_t) * n);
  sfft_output mm;
  GenerateModeMap(n, k, &mm);
  GenerateXhatAlt(n, mm, sigma, xh);

  fftw_plan f_plan = fftw_plan_dft_1d(n, xh, x, FFTW_BACKWARD, FFTW_ESTIMATE);
  fftw_execute_dft(f_plan, xh, x);

  sfft_output out;
  sfft_plan *s_plan = sfft_make_plan(n, k, version, FFTW_ESTIMATE);
  while (state.KeepRunning()) {
    ::sfft_exec(s_plan, x, &out);
  }

  fftw_destroy_plan(f_plan);
  sfft_free(x);
  sfft_free(xh);
}

static void CustomArguments(benchmark::internal::Benchmark *b) {
  const int list_ver[] = {SFFT_VERSION_1, SFFT_VERSION_2};
  for (int j = 0; j < 2; ++j) {
    const int ver = list_ver[j];
    for (int i = 13; i <= 24; ++i) {
      std::vector<int> inputs;
      inputs.push_back(1 << i);
      inputs.push_back(50);
      inputs.push_back(ver);
      b->Args(inputs);
    }
  }

  const int list_k[] = {50, 100, 200, 500, 1000, 2000, 2500, 4000};
  for (int j = 0; j < 2; ++j) {
    const int ver = list_ver[j];
    for (int i = 0; i < 8; ++i) {
      std::vector<int> inputs;
      inputs.push_back(1 << 22);
      inputs.push_back(list_k[i]);
      inputs.push_back(ver);
      b->Args(inputs);
    }
  }
}
BENCHMARK(BM_SFFT)->Apply(CustomArguments);

} // namespace mps

int main(int argc, char **argv) {
  google::InitGoogleLogging(argv[0]);
  google::InstallFailureSignalHandler();
  benchmark::Initialize(&argc, argv);
  benchmark::RunSpecifiedBenchmarks();
}