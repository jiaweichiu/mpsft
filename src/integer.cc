/*
 * Copyright (c) 2017 Jiawei Chiu
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301 USA
 *
 */
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