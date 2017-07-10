# MPSFT: Matrix Pencil Sparse Fourier Transform

For now, we will use the GPL license:
https://www.gnu.org/licenses/gpl-3.0.en.html

This is mainly due to our use of FFTW. We do intend to move away from that in the near future, say using KissFFT.

If you use our results, please reference
https://dspace.mit.edu/handle/1721.1/83691?show=full

It has been a long time. I need to update some details, and at least put on arxiv. For now, the above is all we have.

## Install FFTW

```
./configure --enable-shared --enable-threads --enable-openmp
make
make install
```