HOME=/home/jchiu
BOOSTBASE=/usr/local
FFTWBASE=/usr/local
#EIGENBASE=$(HOME)/eigen-eigen-5097c01bcdc4

#HOME=/home/jiawei
#BOOSTBASE=/usr/local
#FFTWBASE=/usr/local
#EIGENBASE=$(HOME)/eigen-eigen-5097c01bcdc4

#LBOOST=$(BOOSTBASE)/lib/libboost_system.a $(BOOSTBASE)/lib/libboost_timer.a $(BOOSTBASE)/lib/libboost_chrono.a -lrt
LBOOST=-lboost_system -lboost_filesystem -lboost_timer -lboost_chrono
# LFFTW=$(FFTWBASE)/lib/libfftw3.a
LFFTW=-lfftw3

# HEADER=-I$(EIGENBASE) -I$(BOOSTBASE)/include -I$(FFTWBASE)/include
HEADER=-I$(BOOSTBASE)/include -I$(FFTWBASE)/include
LIB=$(LBOOST) $(LFFTW) -lglog

CC=g++
RM=rm -f
#CXXFLAGS=-Wall -Wconversion -g
CXXFLAGS=-Wall -Wconversion -O3 -mtune=native -ffast-math -fopenmp-simd $(HEADER)
LDFLAGS=-g
#LDLIBS=-lfftw3 -lfftw3f -lm -lrt
LDLIBS=-lm -lrt  $(LIB) 

DESTS=experiment generate_graphs kaiserbessel
SRCS:=$(filter-out $(DESTS:%=%.cc), $(wildcard *.cc))
OBJS:=$(subst .cc,.o,$(SRCS))

all: $(DESTS)

$(DESTS): $(OBJS) 

clean:
	$(RM) *.o $(DESTS)
