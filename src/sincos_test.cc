#include "sincos.h"

namespace mps {

TEST_CASE("SinTwoPiBasic", "") {
  const Int n = 1234;
  for (Int i = 0; i <= n; ++i) {
    const double x = double(i) / double(n);
    REQUIRE(SinTwoPi(x) == Approx(std::sin(2 * M_PI * x)));
  }
}

} // namespace mps