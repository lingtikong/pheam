#ifndef TIMER_H
#define TIMER_H

#include "stdio.h"
#include "stdlib.h"
#include "time.h"

class Timer {
public:
  Timer();

  void start();
  void stop();
  void print();
  double cpu_time();
  double wall_time();
  double up2now();
  double sincelast();

private:
  clock_t t0, t1, t2;
  double cpu_time_used;

  time_t tbeg, tend;

  int flag;
};

#endif
