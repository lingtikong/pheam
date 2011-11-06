#include "timer.h"

Timer::Timer()
{
  flag = 0;
  start();
return;
}

void Timer::start()
{
 t1 = clock();
 time(&tbeg);

 flag |= 1;

return;
}

void Timer::stop()
{
  if ( flag&1 ) {
    t2 = clock();
    time(&tend);

    flag |= 2;
  }
return;
}

void Timer::print()
{
  if ( (flag&3) != 3) return;

  cpu_time_used = ((double) (t2 - t1)) / CLOCKS_PER_SEC;
  printf("Total CPU time used: %g seconds; walltime: %g seconds.\n", cpu_time_used, difftime(tend,tbeg));

return;
}

double Timer::cpu_time()
{
  if ( (flag&3) != 3) return 0.;
  else return ((double) (t2 - t1)) / CLOCKS_PER_SEC;
}

double Timer::wall_time()
{
  if ( (flag&3) != 3) return 0.;
  else return difftime(tend,tbeg);
}
