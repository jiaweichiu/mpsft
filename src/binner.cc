#include <glog/logging.h>

#include "binner.h"

Binner::Binner(const BinnerOpt& opt) {
  CHECK_GT(opt.n, 0);
  CHECK_GT(opt.b, 0);
  CHECK_GT(opt.delta, 0);

  n = opt.n;
  b = opt.b;
  delta = opt.delta;
}