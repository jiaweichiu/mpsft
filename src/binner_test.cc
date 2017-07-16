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
#include "catch.hpp"

#include "base.h"
#include "binner.h"
#include "window.h"

namespace mps {

TEST_CASE("BinInTime", "") {
  const Int n = 1109;
  const Int bins = 1;
  const Int bits = 2;
  const Int q = 100;

  const ModeMap mm = {{500, Cplex(1.0, 0)}};
  const CplexArray x = EvaluateModes(n, mm);

  Window win(n, bins, 1e-6);
  Transform tf(n, 1, 0, 0);

  BinnerFast binner(win, bits);
  CplexMatrix out_time(1 + 2 * bits, bins);
  binner.BinInTime(x, tf, q, &out_time);

  BinnerSimple binner0(win, bits);
  CplexMatrix out_time0(1 + 2 * bits, bins);
  binner0.BinInTime(x, tf, q, &out_time0);

  for (Int i = 0; i < 1 + 2 * bits; ++i) {
    REQUIRE(std::abs(out_time[i][0] - out_time0[i][0]) == Approx(0));
  }
}

void BinnerBasic(int binner_type) {
  const Int n = 1109;
  const Int bins = 5;
  const Int bits = 2;
  const Int q = 106;
  Window win(n, bins, 1e-6);

  // Prepare x_hat and x.
  const ModeMap mm = {{550, Cplex(2.0, 0)}};
  const CplexArray x = EvaluateModes(n, mm);

  // Prepare binning.
  Transform tf(n, 3, 847, 45);
  std::unique_ptr<Binner> binner(Binner::Create(binner_type, win, bits));

  // BinInTime.
  CplexMatrix out_time(1 + 2 * bits, bins);
  binner->BinInTime(x, tf, q, &out_time);

  REQUIRE(RE(out_time[0][1]) == Approx(1.12652));
  REQUIRE(IM(out_time[0][1]) == Approx(-0.108838));
  for (Int i = 0; i < bins; ++i) {
    if (i != 1) {
      REQUIRE(std::abs((out_time[0])[i]) == Approx(0));
    }
  }

  // BinInFreq.
  CplexMatrix out_freq(1 + 2 * bits, bins);
  out_freq.Clear();
  binner->BinInFreq(mm, tf, q, &out_freq); // Subtract.

  REQUIRE(RE(out_freq[0][1]) == Approx(-1.12652));
  REQUIRE(IM(out_freq[0][1]) == Approx(0.108838));
  for (Int i = 0; i < bins; ++i) {
    if (i != 1) {
      REQUIRE(std::abs(out_freq[0][i]) == Approx(0));
    }
  }

  // Compare out_time and out_freq.
  for (Int i = 0; i < 1 + 2 * bits; ++i) {
    for (Int j = 0; j < bins; ++j) {
      REQUIRE(std::abs(out_time[i][j] + out_freq[i][j]) == Approx(0));
    }
  }
}

TEST_CASE("BinnerBasic_Simple", "") { BinnerBasic(kBinnerSimple); }

TEST_CASE("BinnerBasic_Fast", "") { BinnerBasic(kBinnerFast); }

void BinnerBigger(int binner_type) {
  const Int n = kPrimes[20];
  const Int bins = 5;
  const Int bits = 2;
  const Int q = 0xFFFFFF;
  Window win(n, bins, 1e-6);

  // Prepare x_hat and x.
  const ModeMap mm = {{565336, Cplex(2.0, 1.0)}};
  const CplexArray x = EvaluateModes(n, mm);

  // Prepare binning.
  Transform tf(n, 0x3FFFFFFF, 0xEEEEEEEE, 0xDDDDDD);
  std::unique_ptr<Binner> binner(Binner::Create(binner_type, win, bits));

  CplexMatrix out_time(1 + 2 * bits, bins);
  binner->BinInTime(x, tf, q, &out_time);

  CplexMatrix out_freq(1 + 2 * bits, bins);
  out_freq.Clear();
  binner->BinInFreq(mm, tf, q, &out_freq); // Subtract.

  for (Int i = 0; i < bins; ++i) {
    REQUIRE(std::abs(out_time[0][i] + out_freq[0][i]) == Approx(0));
  }
}

TEST_CASE("BinnerBigger_Simple", "") { BinnerBigger(kBinnerSimple); }

TEST_CASE("BinnerBigger_Fast", "") { BinnerBigger(kBinnerFast); }

} // namespace mps