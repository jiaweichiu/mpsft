#pragma once

// using my_uint = uint64_t;
// using my_sint = int64_t;

using my_uint = uint32_t;
using my_sint = int32_t;

struct magics_info {
  my_sint multiplier; // the "magic number" multiplier
  unsigned shift;     // shift for the dividend after multiplying
};

struct magicu_info {
  my_uint multiplier;  // the "magic number" multiplier
  unsigned pre_shift;  // shift for the dividend before multiplying
  unsigned post_shift; // shift for the dividend after multiplying
  int increment; // 0 or 1; if set then increment the numerator, using one of
                 // the two strategies
};

magics_info compute_signed_magic_info(my_sint D);
magicu_info compute_unsigned_magic_info(my_uint D, unsigned num_bits);

/*
   To emit code for n/d, rounding towards zero, use the following sequence:

     m = compute_signed_magic_info(D)
     emit("result = (m.multiplier * n) >> SINT_BITS");
     if d > 0 and m.multiplier < 0: emit("result += n")
     if d < 0 and m.multiplier > 0: emit("result -= n")
     if m.post_shift > 0: emit("result >>= m.shift")
     emit("result += (result < 0)")
*/

// Divide n by d. Associated with d is a precomputed <multiplier, shift>.
// #pragma omp declare simd
inline int32_t ApplyMagic(int32_t n, int32_t multiplier, int shift) {
  // __int128 x = __int128(n) * __int128(multiplier);
  // x >>= 64;
  int64_t x = int64_t(n) * int64_t(multiplier);
  x >>= 32;
  const int32_t y = (multiplier < 0) ? (x + n) : x;
  return int32_t(y >> shift);
}
