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
#include "integer.h"
#include "window.h"

namespace mps {

// tau = q +/- (1 << b) where 0 <= b < bits.
// Num of rows of "out" is 2*bits+1.

class BinInTime {
public:
  BinInTime(const Window &win, int32_t bits);
  virtual void Run(const CplexArray &x, const Transform &tf, int32_t q,
                   CplexMatrix *out) = 0;
  static BinInTime *Create(int binner_type, const Window &win, int32_t bits);

protected:
  const Window &win_;
  int32_t bits_;
  std::unique_ptr<FFTPlan> plan_;
};

class BinInTimeV0 : public BinInTime {
public:
  BinInTimeV0(const Window &win, int32_t bits);
  void Run(const CplexArray &x, const Transform &tf, int32_t q,
           CplexMatrix *out) override;

private:
  CplexArray scratch_;
};

// Make use of symmetry in tau's to half number of sinusoids.
// NOTE: Focus is on correctness not efficiency.
// Code is already made more complicated than V0 due to the trick to halve the
// number of Sinusoids. Do not complicate it further.
class BinInTimeV1 : public BinInTime {
public:
  BinInTimeV1(const Window &win, int32_t bits);
  void Run(const CplexArray &x, const Transform &tf, int32_t q,
           CplexMatrix *out) override;

private:
  CplexArray scratch_; // Size bins.
  CplexArray scratch2_;
};

// Extend from V1. Do vectorization! We use SinCos approximation.
// :binner_test pass, but this is failing :iterate_test.
// BinInTimeV3 shall use boost::SIMD instead for better precision.
class BinInTimeV2 : public BinInTime {
public:
  BinInTimeV2(const Window &win, int32_t bits);
  void Run(const CplexArray &x, const Transform &tf, int32_t q,
           CplexMatrix *out) override;

private:
  CplexArray s1_; // Size p=win.p().
  CplexArray s2_;
  Int32Array idx1_;
  Int32Array idx2_;
  DoubleArray dbl_;
};

// Extend from V1. Do vectorization! We use SinCos approximation.
// :binner_test pass, but this is failing :iterate_test.
// BinInTimeV3 shall use boost::SIMD instead for better precision.
class BinInTimeV3 : public BinInTime {
public:
  BinInTimeV3(const Window &win, int32_t bits);
  void Run(const CplexArray &x, const Transform &tf, int32_t q,
           CplexMatrix *out) override;

private:
  CplexArray s1_; // Size p=win.p().
  CplexArray s2_;
  Int32Array idx1_;
  Int32Array idx2_;
  DoubleArray d1_;
  DoubleArray d2_;
  DoubleArray d3_;
  DoubleArray d4_;
};

class BinInFreq {
public:
  BinInFreq(const Window &win, int32_t bits);
  virtual void Run(const ModeMap &mm, const Transform &tf, int32_t q,
                   CplexMatrix *out) = 0;
  static BinInFreq *Create(int binner_type, const Window &win, int32_t bits);

protected:
  const Window &win_;
  int bits_;
};

class BinInFreqV0 : public BinInFreq {
public:
  BinInFreqV0(const Window &win, int32_t bits);
  void Run(const ModeMap &mm, const Transform &tf, int32_t q,
           CplexMatrix *out) override;
};

class BinInFreqV1 : public BinInFreq {
public:
  BinInFreqV1(const Window &win, int32_t bits);
  void Run(const ModeMap &mm, const Transform &tf, int32_t q,
           CplexMatrix *out) override;
};

} // namespace mps