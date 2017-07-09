#include <cstdint>

#include "catch.hpp"

#include "perm.h"

constexpr int32_t n = 536870909;  // Prime.

TEST_CASE("PermTest", "Basic") {
  Perm perm(n, 10000000, 10000000);
  REQUIRE(perm.Forward(10000000) == 287006024);
  REQUIRE(perm.Backward(287006024) == 10000000);
}
