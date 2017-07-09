#include "catch.hpp"

#include "base.h"
#include "binner.h"
#include "perm.h"

namespace mps {

constexpr Int n = 536870909; // Prime.

TEST_CASE("BinnerBasic", "") {
  BinnerOpt opt;
  opt.n = n;
  opt.bins = 5;
  opt.delta = 1e-6;
  Window win(opt);
}

} // namespace mps