#include "freqid.h"

namespace mps {

using Mat22 = Eigen::Matrix<Cplex, 2, 2>;

bool MatPencil(Cplex u0, Cplex u1, Cplex u2, Real *sigma) {
  Mat22 a;
  a(0, 0) = u0;
  a(0, 1) = u2;
  a(1, 0) = u1;
  a(1, 1) = u0;

  Eigen::JacobiSVD<Mat22> solver(a, Eigen::ComputeFullV);

  const auto &sv = solver.singularValues();
  sigma[0] = sv[0];
  sigma[1] = sv[1];
  // LOG(INFO) << "sigma=" << sv;

  const Mat22 &v = solver.matrixV();
  // [v(0, 0), v(1, 0)] has unit norm. This is the dominant singular vector.
  // [v(0, 1), v(1, 1)] has unit norm.
  // const Cplex w = std::conj(v(0, 0)) * v(1, 0);
  // LOG(INFO) << w;
  // return IM(w) < 0;
  return IM(std::conj(v(0, 0)) * v(1, 0)) < 0;
}

} // namespace mps