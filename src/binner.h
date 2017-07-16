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
#include "window.h"

namespace mps {

constexpr int kBinnerSimple = 1;
constexpr int kBinnerFast = 2;

class Binner {
public:
  Binner(const Window &win, Int bits);

  // tau = q +/- (1 << b) where 0 <= b < bits.
  // Num of rows of "out" is 2*bits+1.
  virtual void BinInTime(const CplexArray &x, const Transform &tf, Int q,
                         CplexMatrix *out) = 0;
  virtual void BinInFreq(const ModeMap &mm, const Transform &tf, Int q,
                         CplexMatrix *out) = 0;

  static Binner *Create(int binner_type, const Window &win, Int bits);

protected:
  const Window &win_;
  Int bits_;
  std::unique_ptr<FFTPlan> plan_;
  std::unique_ptr<CplexArray> scratch_;
};

class BinnerSimple : public Binner {
public:
  BinnerSimple(const Window &win, Int bits);

  void BinInTime(const CplexArray &x, const Transform &tf, Int q,
                 CplexMatrix *out) override;
  void BinInFreq(const ModeMap &mm, const Transform &tf, Int q,
                 CplexMatrix *out) override;
};

// Make use of symmetry in tau's to half number of sinusoids.
class BinnerFast : public Binner {
public:
  BinnerFast(const Window &win, Int bits);

  void BinInTime(const CplexArray &x, const Transform &tf, Int q,
                 CplexMatrix *out) override;
  void BinInFreq(const ModeMap &mm, const Transform &tf, Int q,
                 CplexMatrix *out) override;

protected:
  std::unique_ptr<CplexArray> scratch2_;
};

} // namespace mps