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

#include <benchmark/benchmark.h>

// Parameters (n, k): n is size of x and there are k modes.
// TODO: Allow sigma to be varied.
static void BM_demo1(benchmark::State &state) {

  const Int n = state.range(0);
  const Int num_modes = 100; // Shouldn't matter.
  const double sigma = 0.1;

  ModeMap mm;
  for (Int i = 0; i < num_modes; ++i) {
    mm[PosMod(RandomInt(), n)] = Cplex(RandomNormal(), RandomNormal());
  }
  const CplexArray xh = GenerateXhat(n, mm, sigma);
  CplexArray x(n);
}