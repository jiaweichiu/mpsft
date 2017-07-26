#pragma once

#include <complex.h>
#include <stdint.h>
#include <stdlib.h>
#include <tr1/unordered_map>

typedef double complex complex_t;

typedef std::tr1::unordered_map<int, complex_t, std::tr1::hash<int> >
    sfft_output;

enum sfft_version { SFFT_VERSION_1, SFFT_VERSION_2, SFFT_VERSION_3 };

struct sfft_plan {
  sfft_version version;
  unsigned int n;
  unsigned int k;
  void *data;
};

sfft_plan *sfft_make_plan(int n, int k, sfft_version version,
                          int fftw_optimization);
void sfft_free_plan(sfft_plan * plan);

void *sfft_malloc(size_t s);
void sfft_free(void *);

void sfft_exec(sfft_plan *plan, complex_t *in, sfft_output *out);