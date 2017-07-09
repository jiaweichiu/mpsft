#pragma once

#include "base.h"
#include "perm.h"

namespace mps {

struct BinnerOpt {
  Int n;     // Assumed to be prime.
  Int bins;  // Number of bins.
  Real delta;
};

struct Transform {
  Int a, b, c;
  // y[t] = x[a*t+c] exp(2*pi*i*b*t/n).
  // yh[a*k+b] = xh[k] exp(2*pi*i*c*k/n).
};

struct TauSet {
  // tau's are: q, q +/- s for s in list_s.
  Int q;
  vector<Int> list_s;
};

class Window {
 public:
  Window(const BinnerOpt& opt);

  // Assume 0 <= t <= p2 where p2=(p-1)/2. Assume p is odd.
  inline Real wt(Int t) const { return wt_[t]; }

  inline Int p() const { return p_; }

 private:
  Int n_;
  Int bins_;
  Real delta_;

  Real width_;
  Real sqrt_c_delta_;
  Real sigma_f_;
  Real sigma_t_;
  Int p_;  // Size of support of window in time domain.

  // We precompute the window in time domain.
  // TODO: It might be more efficient not to store this. Instead, for each t, we
  // can compute the window at t and then iterate over tau's.
  vector<Real> wt_;  // Size is (p-1)/2.
};

void BinInTime(const Window& win, const Transform& tf, const TauSet& taus,
               const vector<Cplex>& x);

}  // namespace mps