#include <glog/logging.h>

#define CATCH_CONFIG_RUNNER
// #define CATCH_CONFIG_MAIN
#include "catch.hpp"

#include "base.h"

int main(int argc, char *const argv[]) {
  mps::MainInit(argc, argv);
  return Catch::Session().run(argc, argv);
}