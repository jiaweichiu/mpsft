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
BM_BinInTime/0/22   64254592 ns   65285985 ns         10
BM_BinInTime/1/22   44763790 ns   45316439 ns         16
BM_BinInTime/2/22   38125950 ns   38496863 ns         18
BM_BinInFreq/0/22     213263 ns     215111 ns       3249
BM_BinInFreq/1/22     112746 ns     113640 ns       6134
```

The first parameter selects the binner. The second parameter is `log2(n)`.

Consider `BinInTime`. From V0 to V1, there is a ~1.4X speedup by exploiting symmetry. From V1 to V2, there is a ~1.2X speedup by using smaller loops.

Consider `BinInFreq`, there is a  ~1.9X speedup by exploiting symmetry.

## Benchmarks for FFTW

As expected, FFTW is a lot faster when `n` is a power of 2. Here we see that it is about 3-4X faster. We will need to compare against those powers of 2.

```
-------------------------------------------------------
Benchmark                Time           CPU Iterations
-------------------------------------------------------
BM_FFTW/512           2811 ns       2810 ns     248972
BM_FFTW/1024          7071 ns       7070 ns      98097
BM_FFTW/2048         16153 ns      16150 ns      43235
BM_FFTW/4096         44135 ns      44124 ns      15861
BM_FFTW/8192        103129 ns     103105 ns       6773
BM_FFTW/16384       225994 ns     225949 ns       3040
BM_FFTW/32768       562571 ns     562435 ns       1245
BM_FFTW/65536      1074076 ns    1073675 ns        638
BM_FFTW/131072     2535639 ns    2535202 ns        279
BM_FFTW/262144     6570041 ns    6567777 ns        105
BM_FFTW/524288    13495382 ns   13493038 ns         54
BM_FFTW/1048576   45847787 ns   45832947 ns         16
BM_FFTW/2097152  109663756 ns  109647774 ns          6
BM_FFTW/4194304  238162189 ns  238126294 ns          3
BM_FFTW/509          18251 ns      18245 ns      38231
BM_FFTW/1021         40767 ns      40758 ns      17168
BM_FFTW/2053         99805 ns      99769 ns       6970
BM_FFTW/4099        216602 ns     216570 ns       3233
BM_FFTW/8191        517183 ns     517107 ns       1319
BM_FFTW/16381      1219395 ns    1218845 ns        571
BM_FFTW/32771      2297161 ns    2296780 ns        305
BM_FFTW/65537      3668541 ns    3667993 ns        191
BM_FFTW/131071    14585802 ns   14583150 ns         48
BM_FFTW/262147    23815421 ns   23812022 ns         30
BM_FFTW/524287    81973846 ns   81962204 ns          8
BM_FFTW/1048573  185319479 ns  185291332 ns          4
BM_FFTW/2097143  384113082 ns  384051338 ns          2
BM_FFTW/4194301  814812411 ns  814561923 ns          1
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

Here are the results. We see that the sparsity has to be around 512 in order for MPSFT to be faster than FFTW.


```
-------------------------------------------------------------
Benchmark                      Time           CPU Iterations
-------------------------------------------------------------
BM_Demo1/4194301/64     85108215 ns   84031097 ns          8
BM_Demo1/4194301/128   101302040 ns  100489462 ns          7
BM_Demo1/4194301/256   146711994 ns  146038058 ns          6
BM_Demo1/4194301/512   222997975 ns  222333162 ns          3
BM_Demo1/4194301/1024  385268980 ns  384430733 ns          2
BM_Demo1/4194301/2048  634506174 ns  633584324 ns          1
```

A couple of parameters are not the most aggressive. For example, `window_delta` can be slightly bigger.

We use `trials=1` because there is quite little noise here and there is no need to do any probability amplification. We expect `trials` to be odd and the running time is roughly proportional to `trials`. So using `trials=3` will slow us down by ~3X.

# Profiling

Install gperftools. Link binary to this.

```shell
bazel build -c opt --config=opt --linkopt="-lprofiler" :demo1_main

BIN=./bazel-bin/demo1_main

$BIN
pprof --svg $BIN /tmp/demo1_main.prof > profile/demo1_main_20170717.svg
pprof --pdf $BIN /tmp/demo1_main.prof > profile/demo1_main_20170717.pdf
pprof --text $BIN /tmp/demo1_main.prof > profile/demo1_main_20170717.txt
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

# Miscel

## CPU scaling

Switch between `performance` and `powersave`.

```shell
cat /sys/devices/system/cpu/cpu*/cpufreq/scaling_governor

for CPUFREQ in /sys/devices/system/cpu/cpu*/cpufreq/scaling_governor; do [ -f $CPUFREQ ] || continue; echo -n performance > $CPUFREQ; done
```