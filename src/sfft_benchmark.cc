/*
bazel build --config=opt2 :sfft_benchmark
*/
#include <benchmark/benchmark.h>
#include <vector>

#include "base.h"
#include "gen.h"
#include "sfft.h"

namespace mps {

static void BM_SFFT(benchmark::State &state) {
  const int32_t n = state.range(0);
  const int32_t k = state.range(1); // Sparsity.
  const sfft_version version = static_cast<sfft_version>(state.range(2));
  constexpr double sigma = 1e-1;

  CplexArray x(n);
  CplexArray xh(n);
  ModeMap mm;

  GenerateModeMap(n, k, &mm);
  // GenerateXhat(n, mm_, sigma, &xh_);
  GenerateXhatAlt(n, mm, sigma, &xh);
  FFTPlan plan(n, FFTW_BACKWARD);
  plan.Run(xh, &x);


  ::sfft_plan *s_plan = sfft_make_plan(n, k, version, FFTW_ESTIMATE);
  ::sfft_output out;
  while (state.KeepRunning()) {
    ::sfft_exec(s_plan, x.data(), &out);
  }
}

static void CustomArguments(benchmark::internal::Benchmark *b) {
  vector<int> list_n = {50, 100, 200, 500, 1000, 2000, 2500, 4000};
  vector<int> list_ver = {SFFT_VERSION_1, SFFT_VERSION_2};
  for (int ver : list_ver) {
    for (int i = 13; i <= 24; ++i) {
      b->Args({1 << i, 50, ver});
    }
    for (int n : list_n) {
      b->Args({1 << 2, n, ver});
    }
  }
}
BENCHMARK(BM_SFFT)->Apply(CustomArguments);

} // namespace mps