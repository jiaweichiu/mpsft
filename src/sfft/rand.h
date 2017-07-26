#pragma once

#include <stdint.h>

namespace mps {

void RandomSeed(int32_t seed);
int32_t RandomInt32();
int64_t RandomInt64();
double RandomNormal();

}  // namespace mps