#include "perm.h"
#include "base.h"

namespace mps {

namespace {

// PowMod returns mod(b^e, m).
// b=base, e=exponent, m=modulus.
Int PowMod(Int b, Int e, Int m) {
  Int r = 1;
  while (e > 0) {
    if (e & 1) { // Odd exponent
      r = Mod(Long(r) * Long(b), m);
    }
    e >>= 1;
    b = Mod(Long(b) * Long(b), m);
  }
  return r;
}

// Returns a's inverse modulo m. Caution: Assumes prime m.
Int InvMod(Int a, Int m) { return PowMod(a, m - 2, m); }

} // namespace

Perm::Perm(Int n, Int a, Int b) : n_(n), a_(a), b_(b) { a_inv_ = InvMod(a, n); }

Perm::Perm(Int n) : n_(n) {
  a_ = (RandomInt() % (n_ - 1)) + 1;
  b_ = RandomInt() % n;
}

} // namespace mps