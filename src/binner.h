#pragma once

#include "perm.h"

struct BinnerOpt {
  int n; // Assumed to be prime.
  int b; // Number of bins.
  double delta;
};

class Binner {
public:
  Binner(const BinnerOpt &opt);

private:
  int n_;
  int b_;
  double delta_;

  double width_;
  double sqrt_c_delta_;
  double sigma_f_;
  double sigma_t_;
};