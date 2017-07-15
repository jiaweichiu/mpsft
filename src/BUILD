cc_library(
  name = "catch",
  srcs = ["catch.cc"],
  hdrs = ["catch.hpp"],
  deps = [
    ":base",
  ],
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
  name = "freqid",
  srcs = ["freqid.cc"],
  hdrs = ["freqid.h"],
  deps = [
    ":base",
  ],
)

cc_test(
  name = "freqid_test",
  srcs = ["freqid_test.cc"],
  deps = [
    ":freqid",
    ":catch",
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
    ":binner",
    ":catch",
  ],
)

cc_library(
  name = "iterate",
  srcs = ["iterate.cc"],
  hdrs = ["iterate.h"],
  deps = [
    ":binner",
    ":freqid",
    ":window",
  ],
)

cc_test(
  name = "iterate_test",
  srcs = ["iterate_test.cc"],
  deps = [
    ":catch",
    ":iterate",
  ],
)
