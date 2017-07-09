#include "perm.h"
#include "rand.h"

namespace {

// PowMod returns mod(b^e, m).
// b=base, e=exponent, m=modulus.
int32_t PowMod(int32_t b, int32_t e, int32_t m) {
  int32_t r = 1;
  while (e > 0) {
    if (e & 1) { // Odd exponent
      r = int32_t((int64_t(r) * int64_t(b)) % m);
    }
    e >>= 1;
    b = int32_t((int64_t(b) * int64_t(b)) % m);
  }
  return r;
}

// Returns a's inverse modulo m. Caution: Assumes prime m.
int32_t InvMod(int32_t a, int32_t m) { return PowMod(a, m - 2, m); }
}

Perm::Perm(int32_t n, int32_t a, int32_t b) : n_(n), a_(a), b_(b) {
  a_inv_ = InvMod(a, n);
}

Perm::Perm(int32_t n) : n_(n) {
  a_ = (RandomInt() % (n_ - 1)) + 1;
  b_ = RandomInt() % n;
}