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
#pragma once

#include <cstdint>

#include "base.h"

namespace mps {

/*
Let's try to work with int32 most of the time, but provide some convenience
functions to avoid overflows.
*/

// Unsafe in the sense that arguments are int32_t and if you put products
// inside, it might overflow without you knowing.
inline int32_t UnsafePosMod1(int32_t x, int32_t d) { return ((x % d) + d) % d; }

inline int32_t UnsafePosMod2(int32_t x, int32_t d) {
  x %= d;
  return x >= 0 ? x : x + d;
}

// Pick the one that is faster.
inline int32_t UnsafePosMod(int32_t x, int32_t d) {
  return UnsafePosMod2(x, d);
}

// Return a * b modulo divisor. Save us from doing a lot of typecasts.
inline int32_t MulMod(int64_t a, int64_t b, int64_t divisor) {
  return (a * b) % divisor;
}

// Same as MulMod but returns positive result.
inline int32_t MulPosMod(int64_t a, int64_t b, int64_t divisor) {
  const int32_t x = MulMod(a, b, divisor);
  return x >= 0 ? x : x + divisor;
}

// Return a * b + c modulo divisor.
inline int32_t MulAddMod(int64_t a, int64_t b, int64_t c, int64_t divisor) {
  return (a * b + c) % divisor;
}

// Same as MulAddMod but returns positive result.
inline int32_t MulAddPosMod(int64_t a, int64_t b, int64_t c, int64_t divisor) {
  const int32_t x = MulAddMod(a, b, c, divisor);
  return x >= 0 ? x : x + divisor;
}

// y[t] = x[a*t+c] exp(2*pi*i*b*t/n).
// yh[a*k+b] = xh[k] exp(2*pi*i*c*k/n).
// Mode permutation: a, b
// Mode modulation: c
// Forward: a*k+b.
// Backward: a_inv*(k-b).
// Very lightweight class. Keep logic in the other code.
struct Transform {
  Transform(int32_t n);
  Transform(int32_t n, int32_t a, int32_t b, int32_t c);
  int32_t a;
  int32_t b;
  int32_t c;
  int32_t a_inv;
};

} // namespace mps