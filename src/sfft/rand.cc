// Do not include these Boost headers with complex.h.
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int_distribution.hpp>
#include <boost/random/normal_distribution.hpp>


boost::random::mt19937 rng;
boost::random::uniform_int_distribution<int32_t> uid32;
boost::random::uniform_int_distribution<int64_t> uid64;
boost::random::normal_distribution<double> nd;

void RandomSeed(int32_t seed) { rng.seed(seed); }
int32_t RandomInt32() { return uid32(rng); }
int64_t RandomInt64() { return uid64(rng); }
double RandomNormal() { return nd(rng); }
