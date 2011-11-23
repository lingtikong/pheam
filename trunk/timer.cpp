#include "timer.h"

Timer::Timer()
{
  flag = 0;
  start();
return;
}

void Timer::start()
{
 t0 = clock(); t1 = t0;
 time(&tbeg);

 flag |= 1;

return;
}

void Timer::stop()
{
  if ( flag&1 ) {
    t2 = clock(); t1 = t2;
    time(&tend);

    flag |= 2;
  }
return;
}

void Timer::print()
{
  if ( (flag&3) != 3) return;

  cpu_time_used = ((double) (t2 - t0)) / CLOCKS_PER_SEC;
  printf("Total CPU time used: %g seconds; walltime: %g seconds.\n", cpu_time_used, difftime(tend,tbeg));

return;
}

double Timer::cpu_time()
{
  if ( (flag&3) != 3) return 0.;
  else return ((double) (t2 - t0)) / CLOCKS_PER_SEC;
}

double Timer::wall_time()
{
  if ( (flag&3) != 3) return 0.;
  else return difftime(tend,tbeg);
}

double Timer::up2now()
{
  if ( (flag&1) != 1) return 0.;
  else {
    t2 = clock(); t1 = t2;
    return ((double) (t2 - t0)) / CLOCKS_PER_SEC;
  }
}

double Timer::sincelast()
{
  if ( (flag&1) != 1) return 0.;
  else {
    t2 = clock();
    return ((double) (t2 - t1)) / CLOCKS_PER_SEC;
    t1 = t2;
  }
}
