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

// tau = q +/- (1 << b) where 0 <= b < bits.
// Num of rows of "out" is 2*bits+1.

class BinInTime {
public:
  BinInTime(const Window &win, Int bits);
  virtual void Run(const CplexArray &x, const Transform &tf, Int q,
                   CplexMatrix *out) = 0;
  static BinInTime *Create(int binner_type, const Window &win, Int bits);

protected:
  const Window &win_;
  Int bits_;
  std::unique_ptr<FFTPlan> plan_;
};

class BinInTimeV0 : public BinInTime {
public:
  BinInTimeV0(const Window &win, Int bits);
  void Run(const CplexArray &x, const Transform &tf, Int q,
           CplexMatrix *out) override;

private:
  CplexArray scratch_;
};

// Make use of symmetry in tau's to half number of sinusoids.
class BinInTimeV1 : public BinInTime {
public:
  BinInTimeV1(const Window &win, Int bits);
  void Run(const CplexArray &x, const Transform &tf, Int q,
           CplexMatrix *out) override;

private:
  CplexArray scratch_; // Size bins.
  CplexArray scratch2_;
};

// Make use of symmetry in tau's to half number of sinusoids.
// Split big loop into multiple small loops.
class BinInTimeV2 : public BinInTime {
public:
  BinInTimeV2(const Window &win, Int bits);
  void Run(const CplexArray &x, const Transform &tf, Int q,
           CplexMatrix *out) override;

private:
  CplexArray scratch_; // Size p=win.p().
  CplexArray scratch2_;
  vector<Int> idx_;
  vector<Int> idx2_;
};

class BinInFreq {
public:
  BinInFreq(const Window &win, Int bits);
  virtual void Run(const ModeMap &mm, const Transform &tf, Int q,
                   CplexMatrix *out) = 0;
  static BinInFreq *Create(int binner_type, const Window &win, Int bits);

protected:
  const Window &win_;
  Int bits_;
};

class BinInFreqV0 : public BinInFreq {
public:
  BinInFreqV0(const Window &win, Int bits);
  void Run(const ModeMap &mm, const Transform &tf, Int q,
           CplexMatrix *out) override;
};

class BinInFreqV1 : public BinInFreq {
public:
  BinInFreqV1(const Window &win, Int bits);
  void Run(const ModeMap &mm, const Transform &tf, Int q,
           CplexMatrix *out) override;
};

} // namespace mps