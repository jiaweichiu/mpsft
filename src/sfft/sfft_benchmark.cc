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
#include <vector>

#include "gen.h"
#include "rand.h"
#include "sfft.h"

namespace mps {

static void BM_SFFT(benchmark::State &state) {
  const int32_t n = state.range(0);
  const int32_t k = state.range(1); // Sparsity.
  const sfft_version version = static_cast<sfft_version>(state.range(2));
  // const double sigma = 1e-8;
  const double sigma = 0;

  complex_t *tmp1 = (complex_t *)sfft_malloc(sizeof(complex_t) * n);
  complex_t *tmp2 = (complex_t *)sfft_malloc(sizeof(complex_t) * n);

  complex_t *x = (complex_t *)sfft_malloc(sizeof(complex_t) * n);
  complex_t *xh = (complex_t *)sfft_malloc(sizeof(complex_t) * n);
  sfft_output mm;
  GenerateModeMap(n, k, &mm);
  CHECK_EQ(mm.size(), k);
  GenerateXhatAlt(n, mm, sigma, xh);

  fftw_plan f_plan = fftw_plan_dft_1d(n, tmp1, tmp2, FFTW_BACKWARD, FFTW_ESTIMATE);
  fftw_execute_dft(f_plan, xh, x);
  
  while (state.KeepRunning()) {
    sfft_plan *s_plan = sfft_make_plan(n, k, version, FFTW_ESTIMATE);
    sfft_output out;
    sfft_exec(s_plan, x, &out);
    sfft_free_plan(s_plan);
  }

  fftw_destroy_plan(f_plan);
  sfft_free(x);
  sfft_free(xh);
  sfft_free(tmp1);
  sfft_free(tmp2);
}

static void CustomArguments(benchmark::internal::Benchmark *b) {
  const int list_ver[] = {SFFT_VERSION_1, SFFT_VERSION_2};
  
  /*for (int j = 0; j < 2; ++j) {
    const int ver = list_ver[j];
    for (int i = 13; i <= 24; ++i) {
      std::vector<int> inputs;
      inputs.push_back(1 << i);
      inputs.push_back(50);
      inputs.push_back(ver);
      b->Args(inputs);
    }
  }*/

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