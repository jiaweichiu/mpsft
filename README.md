# MPSFT: Matrix Pencil Sparse Fourier Transform

## Licensing

For now, we will use the GPL license. This is mainly due to our use of FFTW. We do intend to move away from that in the near future, say using KissFFT.

We are working on a brief paper to put on arXiv. For now, if you use our results, please kindly reference https://dspace.mit.edu/handle/1721.1/83691?show=full

## Installation

### Bazel

We use Bazel (aka Blaze at Google). Do take a look at `.bazelrc` and modify accordingly.

### FFTW

MPSFT will need FFTW to perform FFT on much smaller vectors.

For comparison purposes, it is also good to build FFTW from source optimized for your machine.

```shell
./configure --enable-shared --enable-threads --enable-openmp
make
make install
```

### Eigen

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

You might want to disable CPU scaling. See for example https://askubuntu.com/questions/523640/how-i-can-disable-cpu-frequency-scaling-and-set-the-system-to-performance

## Running

### Benchmarks for binning

```shell
bazel build -c opt --config=opt :binner_bench
./bazel-bin/binner_bench
```