#include "rand.h"

#include <random>

namespace {

std::mt19937 rng;
std::uniform_int_distribution<int32_t> uid;
}

int32_t RandomInt() {
  return uid(rng);
}
