#pragma once

struct BinnerOpt {
  int n;  // Assumed to be prime.
  int b;  // Number of bins.
  double delta;
};

class Binner {
 public:
  Binner(const BinnerOpt& opt);

 private:
  int n;
  int b;
  double delta;

  double c_delta;
};