#pragma once

#include <cstdint>

class Perm {
public:
  Perm(int32_t n, int32_t a, int32_t b);
  Perm(int32_t n);  // Random permutation.

  inline int32_t Forward(int32_t x) const {
    return (int64_t(a_) * int64_t(x) + int64_t(b_)) % n_;
  }

  inline int32_t PosForward(int32_t x) const { return (Forward(x) + n_) % n_; }

  inline int32_t Backward(int32_t x) const {
    return (int64_t(a_inv_) * int64_t(x - b_)) % n_;
  }

  inline int32_t PosBackward(int32_t x) const {
    return (Backward(x) + n_) % n_;
  }

  inline int32_t n() const { return n_; }
  inline int32_t a() const { return a_; }
  inline int32_t b() const { return b_; }

  // Random permutation: k -> ak+b.
  // Random modulation exp(2*pi*i*ck/n)
  int32_t n_;
  int32_t a_;
  int32_t b_;
  int32_t a_inv_;
};