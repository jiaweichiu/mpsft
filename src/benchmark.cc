#include <benchmark/benchmark.h>

#include "base.h"

int main(int argc, char **argv) {
  mps::MainInit(argc, argv);
  benchmark::Initialize(&argc, argv);
  benchmark::RunSpecifiedBenchmarks();
}