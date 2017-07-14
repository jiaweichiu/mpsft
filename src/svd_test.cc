#include "catch.hpp"

#include "base.h"
#include "svd.h"

namespace mps {

TEST_CASE("SVDBasic", "") {
  Mat22 a;
  a(0, 0) = Cplex(-2.3, 5.2);
  a(0, 1) = Cplex(4.0, 3.0);
  a(1, 0) = Cplex(-5.5, -3.1);
  a(1, 1) = a(0, 0);

  Real sigma[2];
  bool r = SVD22(a, sigma);
  REQUIRE(sigma[0] == Approx(11.34419934));
  REQUIRE(sigma[1] == Approx(0.91057194));
  REQUIRE_FALSE(r);
}

TEST_CASE("SVDBasic2", "") {
  Mat22 a;
  a(0, 0) = Cplex(2.3, 5.2);
  a(0, 1) = Cplex(4.0, 3.0);
  a(1, 0) = Cplex(5.5, -3.1);
  a(1, 1) = a(0, 0);

  Real sigma[2];
  bool r = SVD22(a, sigma);
  REQUIRE(sigma[0] == Approx(9.80673965));
  REQUIRE(sigma[1] == Approx(5.77476038));
  REQUIRE(r);
}

TEST_CASE("SVDBasic3", "") {
  Mat22 a;
  a(0, 0) = Cplex(-15.3, -5.2);
  a(0, 1) = Cplex(-1.5, 3.5);
  a(1, 0) = Cplex(8.0, 3.0);
  a(1, 1) = a(0, 0);

  Real sigma[2];
  bool r = SVD22(a, sigma);
  REQUIRE(sigma[0] == Approx(21.27794592));
  REQUIRE(sigma[1] == Approx(12.53032391));
  REQUIRE_FALSE(r);
}

} // namespace mps