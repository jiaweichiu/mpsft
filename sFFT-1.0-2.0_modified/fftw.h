#ifndef FFTW_H
#define FFTW_H

#include "fft.h"
#include <fftw3.h>

#ifdef USE_FLOAT

#define fftw_plan_dft_1d fftwf_plan_dft_1d
#define fftw_plan fftwf_plan
#define fftw_execute fftwf_execute
#define fftw_destroy_plan fftwf_destroy_plan
#define fftw_free fftwf_free

#endif

#define OPTIMIZE_FFTW 0
int fftw_dft(complex_t *out, int n, complex_t *x, int backwards = 0);

#endif
