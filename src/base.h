#pragma once

#include <cmath>
#include <complex>
#include <cstdint>
#include <glog/logging.h>
#include <iostream>

namespace mps {

using Int = int32_t;
using Long = int64_t;
using Real = double;
using Cplex = std::complex<Real>;

Int RandomInt();

inline Int Round(Real x) { return ::floor(x + 0.5); }

} // namespace mps