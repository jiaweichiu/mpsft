#pragma once

#include <glog/logging.h>
#include <cmath>
#include <complex>
#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

namespace mps {

using Int = int32_t;
using Long = int64_t;
using Real = double;
using Cplex = std::complex<Real>;

using CplexPair = std::pair<Cplex, Cplex>;

using std::string;
using std::vector;

Int RandomInt();

// Force cast into longs.
inline Long PosMod(Long x, Long n) { return ((x % n) + n) % n; }
inline Long Mod(Long x, Long n) { return x % n; }

inline Int Round(Real x) { return ::floor(x + 0.5); }
inline Cplex Sinusoid(Real theta) { return Cplex(::cos(theta), ::sin(theta)); }

// Equivalent to multiplying by i: (x+iy)*i = -y+ix.
inline Cplex RotateForward(Cplex x) { return Cplex(-std::imag(x), std::real(x)); }

// Multiply by -i.
inline Cplex RotateBackward(Cplex x) { return Cplex(std::imag(x), -std::real(x)); }

}  // namespace mps