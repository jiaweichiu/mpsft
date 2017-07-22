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

/*
Used for generating test cases. No need to be too efficient. No need to
vectorize. Precision is more important.
*/

#include "base.h"

namespace mps {

// Generate ModeMap with k unique modes, each of magnitude one.
void GenerateModeMap(int32_t n, int32_t k, ModeMap *mm);

void EvaluateModes(int32_t n, const ModeMap &mm, CplexArray *out);

// Add ambience noise such that in the *time domain*, each sample point is
// contaminated by N(0, sigma).
// Note: x(t) = sum_k xh[k] exp(2*pi*i*k*t). This is unnormalized.
// If xh[k] ~ N(0, s*s), then x(t) ~ N(0, s*s*n) where s*s*n=sigma*sigma.
void GenerateXhat(int32_t n, const ModeMap &mm, double sigma, CplexArray *out);

} // namespace mps