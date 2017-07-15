#pragma once

#include <eigen3/Eigen/Dense>

#include "base.h"

namespace mps {

// u2 is u[-1], u1 is u[1], u0 is u[0].
// Returns true if we think the angle is in the second half.
// Assume sigma is an array of size >= 2.
bool MatPencil(Cplex u0, Cplex u1, Cplex u2, double *sigma);

} // namespace mps
