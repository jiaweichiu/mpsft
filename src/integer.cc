#include "integer.h"

namespace mps {

Transform::Transform(Int n) {
  a = (RandomInt32() % (n - 1)) + 1;
  b = RandomInt32() % n;
  c = RandomInt32() % n;
  a_inv = InvMod(a, n);
}

Transform::Transform(Int n, Int a, Int b, Int c) {
  this->a = a;
  this->b = b;
  this->c = c;
  this->a_inv = InvMod(a, n);
}

} // namespace mps