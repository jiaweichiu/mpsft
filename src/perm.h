#pragma once

#include "base.h"

namespace mps {

class Perm {
 public:
  Perm(Int n, Int a, Int b);
  Perm(Int n);  // Random permutation.

  inline Int Forward(Int x) const {
    return (Long(a_) * Long(x) + Long(b_)) % n_;
  }

  inline Int PosForward(Int x) const { return (Forward(x) + n_) % n_; }

  inline Int Backward(Int x) const {
    return (Long(a_inv_) * Long(x - b_)) % n_;
  }

  inline Int PosBackward(Int x) const { return (Backward(x) + n_) % n_; }

  inline Int n() const { return n_; }
  inline Int a() const { return a_; }
  inline Int b() const { return b_; }

  // Random permutation: k -> ak+b.
  // Random modulation exp(2*pi*i*ck/n)
  Int n_;
  Int a_;
  Int b_;
  Int a_inv_;
};

}  // namespace mps