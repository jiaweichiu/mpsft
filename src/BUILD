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
  deps = [
    ":magic",
  ],
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
  name = "integer",
  srcs = ["integer.cc"],
  hdrs = ["integer.h"],
  deps = [
    ":base",
  ],
)

cc_library(
  name = "gen",
  srcs = ["gen.cc"],
  hdrs = ["gen.h"],
  deps = [
    ":base",
    ":integer",
    ":sincos",
  ],
)

cc_test(
  name = "gen_test",
  srcs = ["gen_test.cc"],
  deps = [
    ":catch",
    ":base",
    ":gen",
    ":integer",
  ],
)


cc_binary(
  name = "base_bench",
  srcs = ["base_bench.cc"],
  deps = [
    ":base",
    ":benchmark",
  ],
)

cc_library(
  name = "sincos",
  srcs = ["sincos.cc"],
  hdrs = ["sincos.h"],
  deps = [
    ":base",
  ],
)

cc_binary(
  name = "sincos_bench",
  srcs = ["sincos_bench.cc"],
  deps = [
    ":sincos",
    ":benchmark",
  ],
)

cc_binary(
  name = "sincos_test",
  srcs = ["sincos_test.cc"],
  deps = [
    ":sincos",
    ":catch",
  ],
)

cc_binary(
  name = "integer_bench",
  srcs = ["integer_bench.cc"],
  deps = [
    ":base",
    ":integer",
    ":benchmark",
  ],
)

cc_test(
  name = "integer_test",
  srcs = ["integer_test.cc"],
  deps = [
    ":base",
    ":integer",
    ":catch",
  ],
)

cc_binary(
  name = "fftw_bench",
  srcs = ["fftw_bench.cc"],
  deps = [
    ":base",
    ":gen",
    ":benchmark",
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
    ":sincos",
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
    ":integer",
    ":window",
    ":sincos",
  ],
)

cc_test(
  name = "binner_test",
  srcs = ["binner_test.cc"],
  deps = [
    ":binner",
    ":gen",
    ":catch",
  ],
)

cc_binary(
  name = "binner_bench",
  srcs = ["binner_bench.cc"],
  deps = [
    ":binner",
    ":gen",
    ":benchmark",
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

cc_library(
  name = "benchmark_lib",
  srcs = ["libbenchmark.a"],
)

cc_library(
  name = "benchmark",
  srcs = ["benchmark.cc"],
  deps = [
    ":base",
    ":benchmark_lib",
  ],
)

cc_library(
  name = "demo1",
  srcs = ["demo1.cc"],
  hdrs = ["demo1.h"],
  deps = [
    ":base",
    ":iterate",
  ],
)

cc_binary(
  name = "demo1_main",
  srcs = ["demo1_main.cc"],
  deps = [
    ":base",
    ":iterate",
    ":demo1",
  ],
  linkopts = ["-lprofiler"],
)

cc_binary(
  name = "demo1_bench",
  srcs = ["demo1_bench.cc"],
  deps = [
    ":base",
    ":iterate",
    ":demo1",
    ":benchmark",
  ],
)

cc_library(
  name = "magic",
  srcs = ["magic.cc"],
  hdrs = ["magic.h"],
)
