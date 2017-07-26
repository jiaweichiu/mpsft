#ifndef COMPUTEFOURIER_H
#define COMPUTEFOURIER_H

#include "fft.h"

#include "fftw.h"
#include "filters.h"
#include <complex.h>
#include <map>

#define OPTIMIZE_FFTW 0
//#define  WITH_COMB 0

extern bool WITH_COMB;
extern bool ALGORITHM1;
extern bool VERBOSE;
extern bool TIMING;

// Comments located in the cc file.
std::map<int, complex_t>
outer_loop(const complex_t *origx, int n, const Filter &filter,
           const Filter &filter_Est, int B2, int num, int B, int W_Comb,
           int Comb_loops, int loop_threshold, int location_loops, int loops);

#endif
