#pragma once

// Equal-quantile (equal-population) instLumi binning: 10 bins, ~2.64M entries each.
// Derived from rootfiles/Dzero_260426-yrefmva_PbPbUPC_HIForward0_Dpt-2_withlumi.root
// (26,419,060 entries, addlumi.cc output; instLumi in 10^33 cm^-2 s^-1, from instlumi.py).

#include <vector>

namespace params {
  const std::vector<double> lumibins = {
    0.,
    1.5059955e-06,
    1.785839e-06,
    2.105541e-06,
    2.4701735e-06,
    2.8705415e-06,
    3.3701005e-06,
    4.0106025e-06,
    4.896365e-06,
    6.2428255e-06,
    6.581766e-06,
  };
}
