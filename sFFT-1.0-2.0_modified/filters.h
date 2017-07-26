#ifndef GAUSSFILTER_H
#define GAUSSFILTER_H

#include "fft.h"

struct Filter {
  complex_t *time;
  int sizet;
  complex_t *freq;
};

/*
  Create a window function such that:
      the main lobe has width 2 * n * filterfrac
      outside the main lobe, the linf residual is tolerance
  Computes the required w.
  Allocates and returns the filter.
 */
complex_t *make_dolphchebyshev_t(double lobefrac, double tolerance, int &w);

complex_t *make_gaussian_t(double lobefrac, double tolerance, int &w);

complex_t *make_kaiserbessel_t(double lobefrac, double tolerance, int &w);

/*
  Modifies a w-dimensional window function to have n-dimensional FFT
  the sum of b adjacent ones previously.

  Allocates and returns a Filter instance pointing to the modified x and an
  n-dimensional FFT of it.
 */

Filter make_multiple_t(complex_t *x, int w, int n, int b);

#endif
