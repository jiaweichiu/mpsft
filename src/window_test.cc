#include "catch.hpp"

#include "base.h"
#include "window.h"

namespace mps {

constexpr Int n = 1109;
constexpr Int bins = 5;
constexpr Real delta = 1e-4;

TEST_CASE("WindowBasic", "") {
  Window window(n, bins, delta);

  const Int p = window.p();
  const Int p2 = (p - 1) / 2;
  // LOG(INFO) << window.Energy();
}

} // namespace mps