#include <random>

#include "base.h"

namespace mps {

namespace {

std::mt19937 rng;
std::uniform_int_distribution<Int> uid;
}

Int RandomInt() {
  return uid(rng);
}

}  // namespace mps