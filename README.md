<!-- MarkdownTOC -->

- Licensing
- Installation
  - Bazel
  - FFTW
  - Eigen
- Benchmarks
  - Benchmarks for binning
  - Benchmarks for FFTW
  - Benchmarks for MPSFT
- Profiling
- Miscel
  - CPU scaling

<!-- /MarkdownTOC -->

# Licensing

For now, we will use the GPL license. This is mainly due to our use of FFTW. We do intend to move away from that in the near future, say using KissFFT.

We are working on a brief paper to put on arXiv. For now, if you use our results, please kindly reference https://dspace.mit.edu/handle/1721.1/83691?show=full

# Installation

## Bazel

We use Bazel (aka Blaze at Google). Do take a look at `.bazelrc` and modify accordingly.

## FFTW

MPSFT will need FFTW to perform FFT on much smaller vectors.

For comparison purposes, it is also good to build FFTW from source optimized for your machine.

```shell
./configure --enable-shared --enable-threads --enable-openmp
make
make install
```

## Eigen

We use Eigen though this dependency should probably be removed. It would be convenient to `make install` so that the headers are more accessible.

### Install Google benchmarks (optional)

We run benchmarks using Google benchmarks. Remember to build using release mode.

```shell
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ../
make
sudo make install
```

After that, make a symbolic link to the static lib in the `src` directory. You should be able to build `:benchmark` as a quick test.

```shell
ln -s /usr/local/lib/libbenchmark.a libbenchmark.a
```

# Benchmarks

Make sure CPU scaling is turned off. See miscel section. Make sure everything is compiled from source and optimized for your machine. Make sure we use only a single core for fair comparison. (It is easy to parallelize. We do that later.)

## Benchmarks for binning

```shell
bazel build -c opt --config=opt :binner_bench
./bazel-bin/binner_bench
```

```
---------------------------------------------------------
Benchmark                  Time           CPU Iterations
---------------------------------------------------------
BM_BinInTime/0/22   43343983 ns   43346872 ns         16
BM_BinInTime/1/22   31140300 ns   31142036 ns         22
BM_BinInTime/2/22   21985368 ns   21985125 ns         32
BM_BinInFreq/0/22      78472 ns      78477 ns       8893
BM_BinInFreq/1/22      44098 ns      44101 ns      15864
```

The first parameter selects the binner. The second parameter is `log2(n)`.

From V0 to V1 for both `BinInTime` and `BinInFreq`, we exploit symmetry to roughly halve the number of trigonometric operations.

From V1 to V2 and `BinInTime`, we use smaller loops in hope of better compiler optimization, say vectorization, AVX etc.

## Benchmarks for FFTW

As expected, FFTW is a lot faster when `n` is a power of 2. Here we see that it is about 3-4X faster. We will need to compare against those powers of 2.

```
-------------------------------------------------------
Benchmark                Time           CPU Iterations
-------------------------------------------------------
BM_FFTW/512           2808 ns       2807 ns     246425
BM_FFTW/1024          7178 ns       7177 ns      97796
BM_FFTW/2048         16479 ns      16478 ns      42291
BM_FFTW/4096         45325 ns      45323 ns      15458
BM_FFTW/8192        112708 ns     112710 ns       6168
BM_FFTW/16384       230591 ns     230573 ns       3016
BM_FFTW/32768       577172 ns     577180 ns       1206
BM_FFTW/65536      1140694 ns    1140611 ns        601
BM_FFTW/131072     2798055 ns    2797629 ns        253
BM_FFTW/262144     7015315 ns    7015529 ns         96
BM_FFTW/524288    14596828 ns   14596991 ns         49
BM_FFTW/1048576   47893827 ns   47891220 ns         15
BM_FFTW/2097152  116089786 ns  116075660 ns          6
BM_FFTW/4194304  241538333 ns  241533821 ns          3
BM_FFTW/509          18738 ns      18738 ns      37353
BM_FFTW/1021         43045 ns      43046 ns      16243
BM_FFTW/2053        102193 ns     102196 ns       6795
BM_FFTW/4099        224428 ns     224425 ns       3120
BM_FFTW/8191        518752 ns     518746 ns       1322
BM_FFTW/16381      1295849 ns    1295800 ns        537
BM_FFTW/32771      2470733 ns    2470532 ns        285
BM_FFTW/65537      4041268 ns    4041056 ns        156
BM_FFTW/131071    14756750 ns   14756669 ns         47
BM_FFTW/262147    25712783 ns   25711094 ns         27
BM_FFTW/524287    87817229 ns   87809923 ns          8
BM_FFTW/1048573  196519341 ns  196502387 ns          4
BM_FFTW/2097143  406672788 ns  406631881 ns          2
BM_FFTW/4194301  858727800 ns  858697804 ns          1
```

## Benchmarks for MPSFT

We use the following parameters.

```
n = 2^22
window_delta = 1e-5
window_threshold = 0.1
max_stale_iter = 5
min_bins = 101
sigma = 1e-2
```

Here are the results. We see that the sparsity has to be around 1000 in order for MPSFT to be faster than FFTW.

```shell
bazel build -c opt --config=opt :demo1_bench
./bazel-bin/demo1_bench
```

Results:

```
-------------------------------------------------------------
Benchmark                      Time           CPU Iterations
-------------------------------------------------------------
BM_Demo1/4194301/64     48057441 ns   48057968 ns         15
BM_Demo1/4194301/128    61008096 ns   61012363 ns         11
BM_Demo1/4194301/256    82409613 ns   82415425 ns          9
BM_Demo1/4194301/512   128505852 ns  128504101 ns          5
BM_Demo1/4194301/1024  214738172 ns  214702456 ns          3
BM_Demo1/4194301/2048  374692767 ns  374709328 ns          2
```

A couple of parameters are not the most aggressive. For example, `window_delta` can be slightly bigger.

We use `trials=1` because there is quite little noise here and there is no need to do any probability amplification. We expect `trials` to be odd and the running time is roughly proportional to `trials`. So using `trials=3` will slow us down by ~3X.

# Profiling

Install gperftools. Link binary to this.

```shell
bazel build --config=opt --linkopt="-lprofiler" :demo1_main

BIN=./bazel-bin/demo1_main

$BIN
pprof --svg $BIN /tmp/demo1_main.prof > profile/demo1_main_20170718.svg
pprof --pdf $BIN /tmp/demo1_main.prof > profile/demo1_main_20170718.pdf
pprof --text $BIN /tmp/demo1_main.prof > profile/demo1_main_20170718.txt

pprof $BIN /tmp/demo1_main.prof
```

We see that `BinInTime` takes up most of the running time.

```
Total: 164 samples
      81  49.4%  49.4%      152  92.7% mps::BinInTimeV2::Run
      54  32.9%  82.3%       57  34.8% sincos
      14   8.5%  90.9%       14   8.5% __muldc3
       4   2.4%  93.3%        4   2.4% __fmod_finite
       3   1.8%  95.1%        3   1.8% nearbyint
       2   1.2%  96.3%        2   1.2% gammal
       1   0.6%  97.0%        2   1.2% Eigen::internal::svd_precondition_2x2_block_to_be_real::run
       1   0.6%  97.6%        1   0.6% __nss_passwd_lookup
       1   0.6%  98.2%        1   0.6% apply@1ee10
       1   0.6%  98.8%        1   0.6% boost::math::detail::erf_imp
       1   0.6%  99.4%        1   0.6% fmod
       1   0.6% 100.0%        8   4.9% mps::BinInFreqV1::Run
       0   0.0% 100.0%        2   1.2% Eigen::JacobiSVD::compute
       0   0.0% 100.0%      164 100.0% __libc_start_main
```

Here is the [visualization](src/profile/demo1_main_20170717.pdf).

After tweaking sin-cos operations, we have:

```
Total: 96 samples
      71  74.0%  74.0%       89  92.7% mps::BinInTimeV2::Run
      18  18.8%  92.7%       18  18.8% mps::SinCosTwoPi
       2   2.1%  94.8%        2   2.1% Eigen::internal::real_2x2_jacobi_svd
       2   2.1%  96.9%        2   2.1% t2_25
       1   1.0%  97.9%        3   3.1% apply@1b6f0
       1   1.0%  99.0%        4   4.2% mps::BinInFreqV1::Run
       1   1.0% 100.0%        1   1.0% mps::Window::Window
       0   0.0% 100.0%        2   2.1% Eigen::JacobiSVD::compute
```

Here is the [visualization](src/profile/demo1_main_20170718.pdf).

# Miscel

## CPU scaling

Switch between `performance` and `powersave`.

```shell
cat /sys/devices/system/cpu/cpu*/cpufreq/scaling_governor

for CPUFREQ in /sys/devices/system/cpu/cpu*/cpufreq/scaling_governor; do [ -f $CPUFREQ ] || continue; echo -n performance > $CPUFREQ; done
```