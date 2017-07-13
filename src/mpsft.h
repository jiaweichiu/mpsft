#pragma once

#include "base.h"
#include "binner.h"

namespace mps {

void mpsft(const CplexArray &x, Int bins, Real delta, Int trials,
           Real threshold);

} // namespace mps