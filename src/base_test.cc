#include "catch.hpp"

#include "base.h"

namespace mps {

namespace {

void Fill(CplexArray *v) {
  // Just to see if it will crash.
  for (Int i = 0; i < v->Size(); ++i) {
    (*v)[i] = Cplex(0.5, -0.5);
  }
}

} // namespace

TEST_CASE("CplexArrayBasic", "") {
  CplexArray u;
  REQUIRE(u.Size() == 0);
  Fill(&u);

  u.Resize(50);
  REQUIRE(u.Size() == 50);
  Fill(&u);

  CplexArray v(100);
  REQUIRE(v.Size() == 100);
  Fill(&v);

  v.Resize(10);
  REQUIRE(v.Size() == 10);
  Fill(&v);

  v.Resize(10000);
  REQUIRE(v.Size() == 10000);
  Fill(&v);
}

} // namespace mps
