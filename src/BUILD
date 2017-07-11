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

cc_test(
  name = "base_test",
  srcs = ["base_test.cc"],
  deps = [
    ":catch",
    ":base",
  ],
)

cc_library(
  name = "window",
  srcs = ["window.cc"],
  hdrs = ["window.h"],
  deps = [
    ":base",
  ],
)

cc_test(
  name = "window_test",
  srcs = ["window_test.cc"],
  deps = [
    ":catch",
    ":window",
  ],
)


cc_library(
	name = "binner",
	srcs = ["binner.cc"],
	hdrs = ["binner.h"],
  deps = [
    ":base",
    ":window",
  ],
)

cc_test(
  name = "binner_test",
  srcs = ["binner_test.cc"],
  deps = [
    ":catch",
    ":binner",
  ],
)
