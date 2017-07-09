#pragma once

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

}  // namespace mps