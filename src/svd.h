#pragma once

#include <eigen3/Eigen/Dense>

#include "base.h"

namespace mps {

using Mat22 = Eigen::Matrix<Cplex, 2, 2>;

void SVD22(const Mat22& a, vector<Real>* sigma);

} // namespace mps