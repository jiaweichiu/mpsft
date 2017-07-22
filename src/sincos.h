#pragma once

namespace mps {

// NOTE: Approximate result.
// Returns sin(2*pi*x). CAUTION: Assume 0 <= x <= 1.
#pragma omp declare simd
double SinTwoPiApprox(double x);

// NOTE: Approximate result.
// Returns sin(2*pi*x). CAUTION: Assume 0 <= x <= 1.
#pragma omp declare simd
double CosTwoPiApprox(double x);

} // namespace mps