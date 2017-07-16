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
#pragma once

#include <eigen3/Eigen/Dense>

#include "base.h"

namespace mps {

// u2 is u[-1], u1 is u[1], u0 is u[0].
// Returns true if we think the angle is in the second half.
// Assume sigma is an array of size >= 2.
bool MatPencil(Cplex u0, Cplex u1, Cplex u2, double *sigma);

} // namespace mps
