#pragma once

#include <eigen3/Eigen/Dense>

#include "base.h"

namespace mps {

using Mat22 = Eigen::Matrix<Cplex, 2, 2>;

// Returns true if we think the angle is in the second half.
// Assume sigma is an array of size >= 2.
bool SVD22(const Mat22 &a, Real *sigma);

} // namespace mps