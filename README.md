**Table of Contents**  *generated with [DocToc](http://doctoc.herokuapp.com/)*

- [Licensing](#)
- [Installation](#)
  - [Bazel](#)
  - [FFTW](#)
  - [Eigen](#)
    - [Install Google benchmarks (optional)](#)
- [Benchmarks](#)
  - [Benchmarks for binning](#)
  - [Benchmarks for FFTW](#)
  - [Benchmarks for MPSFT](#)
- [Profiling](#)
- [Miscel](#)
  - [CPU scaling](#)

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
BM_BinInTime/1/10     529439 ns     529331 ns       1312
BM_BinInTime/1/13     529241 ns     529011 ns       1317
BM_BinInTime/1/16     526896 ns     526835 ns       1308
BM_BinInTime/2/10     327332 ns     326664 ns       2145
BM_BinInTime/2/13     327918 ns     327265 ns       2142
BM_BinInTime/2/16     331771 ns     331071 ns       2124
BM_BinInFreq/1/10     211107 ns     210678 ns       3312
BM_BinInFreq/1/13     211728 ns     211311 ns       3323
BM_BinInFreq/1/16     210767 ns     210340 ns       3317
BM_BinInFreq/2/10     101143 ns     100943 ns       6913
BM_BinInFreq/2/13     109369 ns     109152 ns       6448
BM_BinInFreq/2/16     110680 ns     110453 ns       6307
```

The first parameter selects the binner. 1 is the basic binner. 2 is the optimized binner. The main optimization has to do with roughly halving the number of trignometric operations.

The second parameter is `log2(n)`.

We see that for `BinInTime`, there is a ~1.6X gain. For `BinInFreq`, there is a ~1.9X gain.

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

Here are the results. We see that the sparsity has to be around 500 or less in order for MPSFT to be faster than FFTW.


```
-------------------------------------------------------------
Benchmark                      Time           CPU Iterations
-------------------------------------------------------------
BM_Demo1/4194301/64     96312215 ns   96303470 ns          7
BM_Demo1/4194301/128   119450439 ns  119445622 ns          6
BM_Demo1/4194301/256   161351308 ns  161341183 ns          4
BM_Demo1/4194301/512   252232483 ns  252212400 ns          3
BM_Demo1/4194301/1024  428662101 ns  428636554 ns          2
BM_Demo1/4194301/2048  750128878 ns  750103387 ns          1
```

A couple of parameters are not the most aggressive. For example, `window_delta` can be slightly bigger.

We use `trials=1` because there is quite little noise here and there is no need to do any probability amplification. We expect `trials` to be odd and the running time is roughly proportional to `trials`. So using `trials=3` will slow us down by ~3X.

# Profiling

Install gperftools. Link binary to this.

```shell
bazel build -c opt --config=opt --linkopt="-lprofiler" :demo1_bench

BIN=./bazel-bin/demo1_bench
CPUPROFILE=/tmp/prof.out $BIN
pprof --web $BIN /tmp/prof.out
```

We see that `BinInTime` takes up most of the running time.

```
Total: 228 samples
      77  33.8%  33.8%      126  55.3% mps::BinnerFast::BinInTime
      42  18.4%  52.2%       43  18.9% sincos
      23  10.1%  62.3%       44  19.3% apply@1e120
      13   5.7%  68.0%       13   5.7% t2_32
      11   4.8%  72.8%       11   4.8% __muldc3
      10   4.4%  77.2%       10   4.4% n1_32
       8   3.5%  80.7%        8   3.5% fftw_cpy2d_pair
       7   3.1%  83.8%        7   3.1% __fmod_finite
       5   2.2%  86.0%        5   2.2% gammal
       4   1.8%  87.7%        4   1.8% fftw_cpy2d
```

# Miscel

## CPU scaling

Switch between `performance` and `powersave`.

```shell
cat /sys/devices/system/cpu/cpu*/cpufreq/scaling_governor

for CPUFREQ in /sys/devices/system/cpu/cpu*/cpufreq/scaling_governor; do [ -f $CPUFREQ ] || continue; echo -n performance > $CPUFREQ; done
```