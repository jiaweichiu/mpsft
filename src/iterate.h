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

#include "base.h"
#include "binner.h"

namespace mps {

struct IterateOptions {
  // bins should scale with number of modes we want to recover, which should
  // decrease exponentially.
  // To reduce error in coefficient estimation, more bins help a lot.
  Int bins = 0;

  // window_delta controls the window. Smaller delta means a more "accurate"
  // window. However, it also means more work for BinInTime.
  double window_delta = 0;

  // Increase trials to improve chance of mode identification.
  // Does not help so much with coefficient estimation.
  Int trials = 0;

  // To speed things up, we do not want to process every bin.
  // If the bin coefficient is too small, we want to skip.
  // Use a larger value to be more conservative and reject more often.
  // If zero, it is ignored.
  double bin_threshold = 0;

  // Each isolated mode is attentuated by some factor. If this factor is too
  // small, we want to skip because coefficient estimation would be more
  // sensitive to errors.
  double window_threshold = 0;

  // Singular value threshold. If the smaller singular value of matrix pencil
  // exceeds this threshold, then there is most likely collision, so skip bin.
  // Use a smaller value to be more conservative and reject more often.
  // If zero, this is ignored.
  double sv_threshold = 0;
};

// x is the signal in time domain.
// mm is the current solution. It is a ModeMap, i.e., map from mode location to
// mode coefficient.
void Iterate(const CplexArray &x, const IterateOptions &opt, ModeMap *mm);

} // namespace mps