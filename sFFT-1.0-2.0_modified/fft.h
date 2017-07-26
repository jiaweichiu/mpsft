#ifndef FFT_H
#define FFT_H

#include <cmath>
#include <complex.h>
#include <vector>

#define Vec(a, b) std::vector<__typeof(*(a))>((a), (a) + (b))

// allow easy change to float or long double
//#define USE_FLOAT
#define USE_DOUBLE

#ifdef USE_FLOAT
typedef float complex complex_t;
typedef float real_t;
#define cexp cexpf
#define exp expf
#endif

#ifdef USE_DOUBLE
typedef double complex complex_t;
typedef double real_t;
#endif

//#define DEBUG

#endif
