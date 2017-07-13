

#include "svd.h"

namespace mps {

void SVD22(const Mat22& a, vector<Real>* sigma) {
  Eigen::JacobiSVD<Mat22> solver(a, Eigen::ComputeFullV);

  const auto& sv = solver.singularValues();
  (*sigma)[0] = sv[0];
  (*sigma)[1] = sv[1];
}

} // namespace mps