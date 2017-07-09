#include <glog/logging.h>

#define CATCH_CONFIG_RUNNER
// #define CATCH_CONFIG_MAIN
#include "catch.hpp"

int main(int argc, char *const argv[]) {
  google::InitGoogleLogging(argv[0]);
  google::InstallFailureSignalHandler();
  return Catch::Session().run(argc, argv);
}