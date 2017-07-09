cc_library(
  name = "catch",
  srcs = ["catch.cc"],
  hdrs = ["catch.hpp"],
)

cc_library(
  name = "rand",
  srcs = ["rand.cc"],
  hdrs = ["rand.h"],
)

cc_library(
  name = "perm",
  srcs = ["perm.cc"],
  hdrs = ["perm.h"],
  deps = [
    ":rand",
  ],
)

cc_binary(
  name = "perm_test",
  srcs = ["perm_test.cc"],
  deps = [
    ":catch",
    ":perm",
  ],
  linkopts = ["-lglog"],
)

cc_library(
	name = "binner",
	srcs = ["binner.cc"],
	hdrs = ["binner.h"],
  deps = [
    ":perm",
  ],
)
