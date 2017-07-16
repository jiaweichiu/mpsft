#include "base.h"

namespace mps {

struct Demo1Options {
  Int trials = 1;
  Int min_bins = 101;
  double window_delta = 1e-6;
  double window_threshold = 0.1;

  // Exit when we haven't done anything for too many iterations.
  Int max_stale_iter = 3;
};

class Demo1 {
public:
  // k is number of modes.
  Demo1(const Demo1Options &opt, Int n, Int k, double sigma);
  void Run();

  // Reports the errors.
  void PostAnalyze();

private:
  const Demo1Options &opt_;
  Int k_;
  CplexArray x_;
  CplexArray xh_;
  double sigma_;

  ModeMap mm_; // Randomly generated.
  ModeMap found_mm_;
};

} // namespace mps