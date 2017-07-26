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

#include "gen.h"
#include "base.h"
#include "integer.h"
#include "sincos.h"

namespace mps {

void GenerateModeMap(int32_t n, int32_t k, ModeMap *mm) {
  mm->clear();
  const size_t kk = k;
  while (mm->size() < kk) {
    const int32_t idx = UnsafePosMod(RandomInt32(), n);
    if (mm->find(idx) != mm->end()) {
      continue;
    }
    Cplex coef(RandomNormal(), RandomNormal());
    coef /= std::abs(coef);
    (*mm)[idx] = coef;
  }
}

void EvaluateModes(int32_t n, const ModeMap &mm, CplexArray *__restrict__ out) {
  CHECK_EQ(out->size(), n);
  out->clear();
  Cplex *data = out->data();

  for (const auto &kv : mm) {
    for (int32_t t = 0; t < n; ++t) {
      const double angle =
          2.0 * M_PI * (double(MulMod(kv.first, t, n)) / double(n));
      data[t] += kv.second * Cplex(::cos(angle), ::sin(angle));
    }
  }
}

// Generate Xhat. Noise energy will sum up to n*sigma*sigma.
void GenerateXhat(int32_t n, const ModeMap &mm, double sigma, CplexArray *out) {
  CHECK_EQ(out->size(), n);
  double noise_energy = 0;
  for (int i = 0; i < n; ++i) {
    (*out)[i] = Cplex(RandomNormal(), RandomNormal());
    noise_energy += AbsSq(RE((*out)[i]), IM((*out)[i]));
  }
  for (const auto &kv : mm) {
    Cplex x = (*out)[kv.first];
    noise_energy -= AbsSq(RE(x), IM(x));
  }
  // Rescale noise by this factor.
  const double factor = sigma / std::sqrt(noise_energy);
  for (int i = 0; i < n; ++i) {
    (*out)[i] *= factor;
  }
  for (const auto &kv : mm) {
    (*out)[kv.first] = kv.second;
  }
}

// Generate Xhat. Noise energy will sum up to n*sigma*sigma.
void GenerateXhatAlt(int32_t n, const ModeMap &mm, double sigma,
                     CplexArray *out) {
  CHECK_EQ(out->size(), n);
  out->clear();
  for (const auto &kv : mm) {
    (*out)[kv.first] = kv.second;
  }
  const double scale = sigma / std::sqrt(2 * n);
  for (int i = 0; i < n; ++i) {
    (*out)[i] += Cplex(RandomNormal(), RandomNormal()) * scale;
  }
}

} // namespace mps