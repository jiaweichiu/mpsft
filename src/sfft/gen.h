#pragma once

#include <stdint.h>

#include "sfft.h"

namespace mps {

void GenerateModeMap(int32_t n, int32_t k, sfft_output *mm);

void GenerateXhatAlt(int32_t n, const sfft_output &mm, double sigma,
                     complex_t *out);

}  // namespace mps