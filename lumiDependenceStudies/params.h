#pragma once

// Equal-quantile (equal-population) instLumi binning: 10 bins, ~2.64M entries each.
// Derived from rootfiles/Dzero_260426-yrefmva_PbPbUPC_HIForward0_Dpt-2_withlumi.root
// (26,419,060 entries, addlumi.cc output; instLumi in 10^33 cm^-2 s^-1, from instlumi.py).

#include <vector>

namespace params {
  // OLD (10 bins, HIForward0 only, leveled events included):
  // const std::vector<double> lumibins = {
  //   0.,
  //   1.5059955e-06,
  //   1.785839e-06,
  //   2.105541e-06,
  //   2.4701735e-06,
  //   2.8705415e-06,
  //   3.3701005e-06,
  //   4.0106025e-06,
  //   4.896365e-06,
  //   6.2428255e-06,
  //   6.581766e-06,
  // };

  // Equal-quantile instLumi binning: 15 bins, ~4.72M entries each (70,782,112 total).
  // Non-leveled events only (lumiLeveled==0, instLumi>0), from the all-PD merged file
  // rootfiles/tuple_Dzero_260426-yrefmva_PbPbUPC_HIForwardMerge_Dpt-2_wLumi.root
  // (instLumi in 10^33 cm^-2 s^-1).
  const std::vector<double> lumibins = {
    0.,
    2.694509930734057e-07,
    1.3488529475580435e-06,
    1.681900471339759e-06,
    1.9474189230095362e-06,
    2.272954588988796e-06,
    2.6008385702880332e-06,
    2.976316409331048e-06,
    3.4008164675469743e-06,
    3.920767539966619e-06,
    4.524284122453537e-06,
    5.112247436045436e-06,
    5.828468601976056e-06,
    6.206209945958108e-06,
    6.292214493441861e-06,
    6.581765774171799e-06,
  };

}
