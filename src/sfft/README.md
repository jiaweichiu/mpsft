CFLAGS="-O3 -DNDEBUG -mtune=native -ffast-math -fopenmp-simd"
LFLAGS="-L/usr/local/lib -lglog -lm -lrt -lpthread -lgomp -lfftw3"

g++ $CFLAGS gen.cc -c

g++ $CFLAGS rand.cc -c

g++ $CFLAGS sfft_benchmark.cc rand.o gen.o \
../libbenchmark.a libsfft.a \
$LFLAGS \
-o sfft_benchmark