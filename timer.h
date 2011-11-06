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

private:
  clock_t t1, t2;
  double cpu_time_used;

  time_t tbeg, tend;

  int flag;
};

#endif
