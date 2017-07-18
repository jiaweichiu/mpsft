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
#include <gperftools/profiler.h>

#include "base.h"
#include "demo1.h"

namespace mps {

void Run(int argc, char *const argv[]) {
  MainInit(argc, argv);
  Demo1Options opt;
  opt.trials = 1;
  opt.min_bins = 201;
  opt.window_delta = 1e-6;
  opt.window_threshold = 0.1;
  opt.max_stale_iter = 5;

  Demo1 demo(opt, kPrimes[24], 1 << 12, 0.1);
  ProfilerStart("/tmp/demo1_main.prof");
  demo.Run();
  ProfilerStop();
  demo.PostAnalyze();
}

} // namespace mps

int main(int argc, char *const argv[]) { mps::Run(argc, argv); }