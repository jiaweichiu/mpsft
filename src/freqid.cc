/*
 * Copyright (c) 2017 Jiawei Chiu
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301 USA
 *
 */
#include "freqid.h"

namespace mps {

using Mat22 = Eigen::Matrix<Cplex, 2, 2>;

bool MatPencil(Cplex u0, Cplex u1, Cplex u2, double *sigma) {
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