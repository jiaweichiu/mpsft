#include <boost/timer/timer.hpp>
#include "MyTimer.h"
#include <stdio.h>

boost::timer::cpu_timer globalTimer;

void Tic() {
  globalTimer.start();
}
double Toc() {
  boost::timer::cpu_times const t(globalTimer.elapsed());
//  boost::timer::nanosecond_type const t2(t.system+t.user);
  boost::timer::nanosecond_type const t2(t.user);
  return double(t2)*1e-9;
}
double TocMore() {
  double t=Toc();
  printf("Time elapsed=%.3e seconds\n",t);
  return t;
}
