#include "filters.h"
#include "plot.h"
#include "utils.h"

#include "fftw.h"
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <cstring>

// Comments located in the header file

double I0(double x) {
  // 1 + sum_L=0^inf  (x/2)^2L  / (L!)^2
  double ans = 1;
  double curval = 1;
  for (int L = 1; curval > 0.001; L++) {
    curval = curval * (x / 2) * (x / 2) / L / L;
    ans += curval;
    // printf("State: %d %lf %lf\n", L, curval, ans);
  }
  return ans;
}

complex_t *make_kaiserbessel_t(double lobefrac, double tolerance, int &w) {
  w = int((1 / M_PI) * (1 / lobefrac) * acosh(1. / tolerance));
  double B = log(1 / tolerance);
  complex_t *x = (complex_t *)malloc(w * sizeof(*x));
  for (int i = 0; i < w; i++) {
    double tmp = (2. * (i - (w - 1) / 2.)) / (w - 1);
    x[i] = I0(B * sqrt(1 - tmp * tmp)) / I0(B);
  }
  return x;
}

double Cheb(double m, double x) {
  if (fabs(x) <= 1)
    return cos(m * acos(x));
  else
    return creal(ccosh(m * cacosh(x)));
}

complex_t *make_dolphchebyshev_t(double lobefrac, double tolerance, int &w) {
  w = int((1 / M_PI) * (1 / lobefrac) * acosh(1. / tolerance));
  if (!(w % 2))
    w--;
  complex_t *x = (complex_t *)malloc(w * sizeof(*x));
  double t0 = cosh(acosh(1 / tolerance) / (w - 1));
  for (int i = 0; i < w; i++) {
    x[i] = Cheb(w - 1, t0 * cos(M_PI * i / w)) * tolerance;
  }
  fftw_dft(x, w, x);
  shift(x, w, w / 2);
  for (int i = 0; i < w; i++)
    x[i] = creal(x[i]);
  return x;
}

complex_t *make_gaussian_t(double lobefrac, double tolerance, int &w) {
  w = int((2 / M_PI) * (1 / lobefrac) * log(1. / tolerance));
  double std_t = (w / 2.) / sqrt(2 * log(1. / tolerance));
  complex_t *x = (complex_t *)malloc(w * sizeof(*x));
  double center = w * 1. / 2;
  for (int i = 0; i < w; i++) {
    real_t dist = fabs(i - center);
    x[i] = exp(-dist * dist / (2 * std_t * std_t));
  }
  return x;
}

double sinc(double x) {
  if (x == 0)
    return 1;
  return sin(x) / (x);
}

Filter make_multiple_t(complex_t *x, int w, int n, int b) {
  assert(b <= n);
  assert(w <= n);
  complex_t *g = (complex_t *)calloc(n, sizeof(*g));
  complex_t *h = (complex_t *)malloc(n * sizeof(*h));
  memcpy(g, x + w / 2, (w - (w / 2)) * sizeof(*g));
  memcpy(g + n - w / 2, x, (w / 2) * sizeof(*g));
  fftw_dft(g, n, g);
  complex_t s = 0;
  for (int i = 0; i < b; i++) {
    s += g[i];
  }
  real_t max = 0;
  int offset = b / 2;
  for (int i = 0; i < n; i++) {
    h[(i + n + offset) % n] = s;
    max = std::max(max, cabs(s));
    s = s + (g[(i + b) % n] - g[i]);
  }
  for (int i = 0; i < n; i++)
    h[i] /= max;

  complex_t offsetc = 1, step = cexp(-2 * M_PI * I * (w / 2) / n);
  for (int i = 0; i < n; i++) {
    // offsetc = cexp(-2*M_PI * I * (w/2) * i / n);
    h[i] *= offsetc;
    offsetc *= step;
  }
  // plot("MOO", map_abs(Vec(h, n)));
  fftw_dft(g, n, h, 1);
  memmove(x, g, w * sizeof(*x));

  /* Figure out exactly the image of the filter?
    memset(g+w, 0, (n-w)*sizeof(*g));
    fftw_dft(h, n, g);
    for(int i = 0; i < n; i++)
    h[i] /= n;
  */

  free(g);
  for (int i = 0; i < w; i++)
    x[i] /= n;

  Filter answer;
  answer.time = x;
  answer.sizet = w;
  answer.freq = h;
  return answer;
}
