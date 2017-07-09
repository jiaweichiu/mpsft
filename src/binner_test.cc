#include "catch.hpp"

#include "base.h"
#include "perm.h"
#include "binner.h"

namespace mps {

constexpr Int n = 536870909; // Prime.

TEST_CASE("BinnerBasic", "") {
  BinnerOpt opt;
  opt.n = n;
  opt.b = 5;
  opt.delta = 1e-6;
  Binner binner(opt);
}

} // namespace mps