#include "fft.h"
#include "fftw.h"
#include "filters.h"
#include "plot.h"
#include "utils.h"
#include <cmath>
#include <cstdlib>
#include <cstring>

int make_box_t(complex_t *x, int w, int n) {
  for (int i = 0; i < w; i++)
    x[i] = 1;
  return 0;
}

int make_hamming_t(complex_t *x, int w, int n) {
  for (int i = 0; i < w; i++)
    x[i] = 0.54 - 0.46 * cos(2 * M_PI * i / (w - 1));
  return 0;
}

int main() {
  int w = 64;
  int w1 = w, w2 = w, w3 = w;
  int n = 8192;
  double frac = 0.05;
  complex_t *filtert = (complex_t *)calloc(n, sizeof(*filtert));
  complex_t *filterf = (complex_t *)calloc(n, sizeof(*filterf));
  complex_t *filtert2 = (complex_t *)calloc(n, sizeof(*filtert));
  complex_t *filterf2 = (complex_t *)calloc(n, sizeof(*filterf));
  complex_t *filtert3 = (complex_t *)calloc(n, sizeof(*filtert));
  complex_t *filterf3 = (complex_t *)calloc(n, sizeof(*filterf));
  complex_t *tmp;

  tmp = make_kaiserbessel_t(frac, 1e-9, w1);
  memcpy(filtert, tmp, w1 * sizeof(*filtert));
  free(tmp);

  tmp = make_dolphchebyshev_t(frac, 1e-9, w3);
  memcpy(filtert3, tmp, w3 * sizeof(*filtert3));
  free(tmp);

  tmp = make_gaussian_t(frac, 1e-9, w2);
  memcpy(filtert2, tmp, w2 * sizeof(*filtert2));
  free(tmp);

  make_multiple_t(filtert, w1, n, 1);
  make_multiple_t(filtert2, w2, n, 1);
  make_multiple_t(filtert3, w3, n, 1);
  fftw_dft(filterf, n, filtert);
  fftw_dft(filterf2, n, filtert2);
  fftw_dft(filterf3, n, filtert3);
  double v1 = cabs(filterf[0]);
  double v2 = cabs(filterf2[0]);
  double v3 = cabs(filterf3[0]);
  v1 /= 1e7;
  v2 /= 1e7;
  v3 /= 1e7;
  for (int i = 0; i < n; i++) {
    filtert[i] /= v1;
    filtert2[i] /= v2;
    filtert3[i] /= v3;
    filterf[i] /= v1;
    filterf2[i] /= v2;
    filterf3[i] /= v3;
  }
  plot("filter time series'\nset logscale y\nz='",
       "Kaiser Bessel\nGaussian\nDolph Chebyshev", map_abs(Vec(filtert, w1)),
       map_abs(Vec(filtert2, w2)), map_abs(Vec(filtert3, w3)));

  plot("filter fourier series'\nset logscale y\nz='",
       "Kaiser Bessel\nGaussian\nDolph Chebyshev", map_abs(Vec(filterf, n)),
       map_abs(Vec(filterf2, n)), map_abs(Vec(filterf3, n)));
}
