#pragma once

#include "base.h"

namespace mps {

struct BinnerOpt {
  Int n; // Assumed to be prime.
  Int b; // Number of bins.
  Real delta;
};

class Binner {
public:
  Binner(const BinnerOpt &opt);

private:
  Int n_;
  Int b_;
  Real delta_;

  Real width_;
  Real sqrt_c_delta_;
  Real sigma_f_;
  Real sigma_t_;
  Int p_;
};

} // namespace mps