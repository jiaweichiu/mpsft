#include <gperftools/profiler.h>

#include "base.h"
#include "demo1.h"

namespace mps {

void Run(int argc, char *const argv[]) {
  MainInit(argc, argv);
  Demo1Options opt;
  opt.trials = 1;
  opt.min_bins = 201;
  opt.window_delta = 1e-6;
  opt.window_threshold = 0.1;
  opt.max_stale_iter = 5;

  Demo1 demo(opt, kPrimes[24], 1 << 12, 0.1);
  ProfilerStart("/tmp/demo1_main.prof");
  demo.Run();
  ProfilerStop();
  demo.PostAnalyze();
}

} // namespace mps

int main(int argc, char *const argv[]) { mps::Run(argc, argv); }