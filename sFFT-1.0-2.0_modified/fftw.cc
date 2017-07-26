#include "fftw.h"
#include <map>

std::map<int, fftw_plan> fftw_plans;

int fftw_dft(complex_t *out, int n, complex_t *x, int backwards) {
  fftw_plan p;
  if (OPTIMIZE_FFTW) { // measure the best plan the first time
    if (fftw_plans.find(n) == fftw_plans.end()) { // no plan made yet
      complex_t *in = (complex_t *)fftw_malloc(sizeof(*in) * n);
      complex_t *out2 = (complex_t *)fftw_malloc(sizeof(*out2) * n);
      p = fftw_plan_dft_1d(
          n, in, out2, backwards ? FFTW_BACKWARD : FFTW_FORWARD, FFTW_MEASURE);
      fftw_plans.insert(std::make_pair(n, p));
      fftw_free(in);
      fftw_free(out2);
    }
  }
  p = fftw_plan_dft_1d(n, x, out, backwards ? FFTW_BACKWARD : FFTW_FORWARD,
                       FFTW_ESTIMATE);
  fftw_execute(p);
  fftw_destroy_plan(p);
  return 0;
}
