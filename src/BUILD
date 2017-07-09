cc_library(
  name = "catch",
  srcs = ["catch.cc"],
  hdrs = ["catch.hpp"],
)

cc_library(
  name = "base",
  srcs = ["base.cc"],
  hdrs = ["base.h"],
)

cc_library(
  name = "perm",
  srcs = ["perm.cc"],
  hdrs = ["perm.h"],
  deps = [
    ":base",
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
