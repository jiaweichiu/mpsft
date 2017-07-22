#include "integer.h"

namespace mps {

namespace {

// PowMod returns mod(b^e, m).
// b=base, e=exponent, m=modulus.
int32_t PowMod(int32_t b, int32_t e, int32_t m) {
  int32_t r = 1;
  while (e > 0) {
    if (e & 1) { // Odd exponent
      r = MulMod(r, b, m);
    }
    e >>= 1;
    b = MulMod(b, b, m);
  }
  return r;
}

// Returns a's inverse modulo m. Caution: Assumes prime m.
int32_t InvMod(int32_t a, int32_t m) { return PowMod(a, m - 2, m); }

} // namespace

Transform::Transform(int32_t n) {
  a = (RandomInt32() % (n - 1)) + 1;
  b = RandomInt32() % n;
  c = RandomInt32() % n;
  a_inv = InvMod(a, n);
}

Transform::Transform(int32_t n, int32_t a, int32_t b, int32_t c) {
  this->a = a;
  this->b = b;
  this->c = c;
  this->a_inv = InvMod(a, n);
}

} // namespace mps