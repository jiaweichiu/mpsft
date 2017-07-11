#pragma once

#include "base.h"

namespace mps {

class Window {
public:
  Window(Int n, Int bins, Real delta);

  // Assume 0 <= t <= p2 where p2=(p-1)/2. Assume p is odd.
  inline Real wt(Int t) const {
    DCHECK_GE(t, 0);
    DCHECK_LE(t, p2());
    return wt_[t];
  }

  inline Int n() const { return n_; }
  inline Int p() const { return p_; }
  inline Int p2() const { return (p_ - 1) / 2; }
  inline Int bins() const { return bins_; }

  // Real Energy() const;

  Real SampleInTime(Int i) const;
  Real SampleInFreq(Real xi) const;

private:
  Int n_;
  Int bins_;
  Real delta_;

  Real width_;
  Real sigma_f_;
  Real sigma_t_;
  Int p_; // Size of support of window in time domain.

  // We precompute the window in time domain.
  // TODO: It might be more efficient not to store this. Instead, for each t, we
  // can compute the window at t and then iterate over tau's.
  vector<Real> wt_; // Size is (p-1)/2.
};

} // namespace mps