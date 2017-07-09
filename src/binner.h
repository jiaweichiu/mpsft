#pragma once

#include "perm.h"

namespace mps {

struct BinnerOpt {
  int32_t n; // Assumed to be prime.
  int32_t b; // Number of bins.
  double delta;
};

class Binner {
public:
  Binner(const BinnerOpt &opt);

private:
  int32_t n_;
  int32_t b_;
  double delta_;

  double width_;
  double sqrt_c_delta_;
  double sigma_f_;
  double sigma_t_;
  int32_t p_;
};

}  // namespace mps