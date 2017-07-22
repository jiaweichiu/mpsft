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
#include "gen.h"
#include "window.h"

namespace mps {

// Check that different versions of BinInTime matches up.
void BinInTimeMatch(int v0, int v1) {
  const int32_t n = 1109;
  const int32_t bins = 1;
  const int32_t bits = 2;
  const int32_t q = 100;
  const ModeMap mm = {{500, Cplex(1.0, 0)}};
  CplexArray x(n);
  EvaluateModes(n, mm, &x);
  Window win(n, bins, 1e-6);
  Transform tf(n, 1, 0, 0);
  std::unique_ptr<BinInTime> binner0(BinInTime::Create(v0, win, bits));
  CplexMatrix out_time0(1 + 2 * bits, bins);
  binner0->Run(x, tf, q, &out_time0);
  std::unique_ptr<BinInTime> binner1(BinInTime::Create(v1, win, bits));
  CplexMatrix out_time1(1 + 2 * bits, bins);
  binner1->Run(x, tf, q, &out_time1);
  for (int32_t i = 0; i < 1 + 2 * bits; ++i) {
    REQUIRE(std::abs(out_time0[i][0] - out_time1[i][0]) == Approx(0));
  }
}
TEST_CASE("BinInTimeMatch_0_1", "") { BinInTimeMatch(0, 1); }
TEST_CASE("BinInTimeMatch_0_2", "") { BinInTimeMatch(0, 2); }

// Check that BinInTime and BinInFreq has identical results.
void TimeFreqMatch(int bin_in_time_type, int bin_in_freq_type) {
  const int32_t n = 1109;
  const int32_t bins = 5;
  const int32_t bits = 2;
  const int32_t q = 106;
  Window win(n, bins, 1e-6);

  // Prepare x_hat and x.
  const ModeMap mm = {{550, Cplex(2.0, 0)}};
  CplexArray x(n);
  EvaluateModes(n, mm, &x);
  Transform tf(n, 3, 847, 45);

  // BinInTime.
  std::unique_ptr<BinInTime> bin_in_time(
      BinInTime::Create(bin_in_time_type, win, bits));
  CplexMatrix out_time(1 + 2 * bits, bins);
  bin_in_time->Run(x, tf, q, &out_time);

  REQUIRE(RE(out_time[0][1]) == Approx(1.12652));
  REQUIRE(IM(out_time[0][1]) == Approx(-0.108838));
  for (int32_t i = 0; i < bins; ++i) {
    if (i != 1) {
      REQUIRE(std::abs((out_time[0])[i]) == Approx(0));
    }
  }

  // BinInFreq.
  std::unique_ptr<BinInFreq> bin_in_freq(
      BinInFreq::Create(bin_in_freq_type, win, bits));
  CplexMatrix out_freq(1 + 2 * bits, bins);
  out_freq.clear();
  bin_in_freq->Run(mm, tf, q, &out_freq); // Subtract.

  REQUIRE(RE(out_freq[0][1]) == Approx(-1.12652));
  REQUIRE(IM(out_freq[0][1]) == Approx(0.108838));
  for (int32_t i = 0; i < bins; ++i) {
    if (i != 1) {
      REQUIRE(std::abs(out_freq[0][i]) == Approx(0));
    }
  }

  // Compare out_time and out_freq.
  for (int32_t i = 0; i < 1 + 2 * bits; ++i) {
    for (int32_t j = 0; j < bins; ++j) {
      REQUIRE(std::abs(out_time[i][j] + out_freq[i][j]) == Approx(0));
    }
  }
}

TEST_CASE("TimeFreqMatch_0_0", "") { TimeFreqMatch(0, 0); }
TEST_CASE("TimeFreqMatch_0_1", "") { TimeFreqMatch(0, 1); }
TEST_CASE("TimeFreqMatch_1_0", "") { TimeFreqMatch(1, 0); }
TEST_CASE("TimeFreqMatch_2_0", "") { TimeFreqMatch(2, 0); }

// Same as TimeFreqMatch but using larger sizes.
void TimeFreqMatchBigger(int bin_in_time_type, int bin_in_freq_type) {
  const int32_t n = kPrimes[20];
  const int32_t bins = 5;
  const int32_t bits = 2;
  const int32_t q = 0xFFFFFF;
  Window win(n, bins, 1e-6);

  // Prepare x_hat and x.
  const ModeMap mm = {{565336, Cplex(2.0, 1.0)}};
  CplexArray x(n);
  EvaluateModes(n, mm, &x);
  Transform tf(n, 0x3FFFFFFF, 0xEEEEEEEE, 0xDDDDDD);

  // BinInTime.
  std::unique_ptr<BinInTime> bin_in_time(
      BinInTime::Create(bin_in_time_type, win, bits));
  CplexMatrix out_time(1 + 2 * bits, bins);
  bin_in_time->Run(x, tf, q, &out_time);

  // BinInFreq.
  std::unique_ptr<BinInFreq> bin_in_freq(
      BinInFreq::Create(bin_in_freq_type, win, bits));
  CplexMatrix out_freq(1 + 2 * bits, bins);
  out_freq.clear();
  bin_in_freq->Run(mm, tf, q, &out_freq); // Subtract.

  // Compare.
  for (int32_t i = 0; i < bins; ++i) {
    REQUIRE(std::abs(out_time[0][i] + out_freq[0][i]) == Approx(0));
  }
}

TEST_CASE("TimeFreqMatchBigger_0_0", "") { TimeFreqMatchBigger(0, 0); }
TEST_CASE("TimeFreqMatchBigger_0_1", "") { TimeFreqMatchBigger(0, 1); }
TEST_CASE("TimeFreqMatchBigger_1_0", "") { TimeFreqMatchBigger(1, 0); }
TEST_CASE("TimeFreqMatchBigger_2_0", "") { TimeFreqMatchBigger(2, 0); }

} // namespace mps