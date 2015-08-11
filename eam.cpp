#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "eam.h"

using namespace std;

#define MAXLINE 512

/*------------------------------------------------------------------------------
 * Constructor is used to
 *----------------------------------------------------------------------------*/
EAM::EAM()
{
  memory = new Memory;
  ntype  = 0;
  initialized = 0;
return;
}

/*------------------------------------------------------------------------------
 * Deconstructor is used to
 *----------------------------------------------------------------------------*/
EAM::~EAM()
{
  memory->destroy(elements);
  memory->destroy(mass);

  delete memory;
}

/*------------------------------------------------------------------------------
 * Method to locate the index of the input element in EAM
 *----------------------------------------------------------------------------*/
int EAM::index(const char *name)
{
  int ip = -1;
  for (int i=0; i<ntype; i++){
    if (strcmp(name, elements[i]) == 0){ ip=i; break;}
  }
  return ip;
}

/*------------------------------------------------------------------------------
 * Method to count # of words in a string, without destroying the string
 *----------------------------------------------------------------------------*/
int EAM::count_words(const char *line)
{
  int n = strlen(line) + 1;
  char *copy;
  memory->create(copy,n,"copy");
  strcpy(copy,line);

  char *ptr;
  if (ptr = strchr(copy,'#')) *ptr = '\0';

  if (strtok(copy," \t\n\r\f") == NULL) {
    memory->destroy(copy);
    return 0;
  }
  n = 1;
  while (strtok(NULL," \t\n\r\f")) n++;

  memory->destroy(copy);
  return n;
}

/*------------------------------------------------------------------------------
 * grab n values from file fp and put them in list.
 * values can be several to one line
 *----------------------------------------------------------------------------*/
void EAM::grab(FILE *fp, int n, double *list)
{
  char *ptr;
  char line[MAXLINE];

  int i = 0;
  while (i < n) {
    fgets(line,MAXLINE,fp);
    ptr = strtok(line," \t\n\r\f");
    list[i++] = atof(ptr);
    while (ptr = strtok(NULL," \t\n\r\f")) list[i++] = atof(ptr);
  }
}

/*------------------------------------------------------------------------------
 * private method to check the equilibrium lattice constant, cohesive energy of
 * BCC and FCC lattice given by the EAM potential(s) read.
 *----------------------------------------------------------------------------*/
void EAM::check_eam()
{
  if (initialized == 0) return;
  printf("\nBased on the EAM potential read:\n");
  for (int i=0; i<44; i++) printf("="); printf("\n");
  for (int i=0; i<15; i++) printf(" "); printf("BCC");
  for (int i=0; i<16; i++) printf(" "); printf("FCC"); printf("\n");
  printf(" Elem  "); for (int i=0; i<17; i++) printf("-");
  printf("   "); for (int i=0; i<17; i++) printf("-"); printf("\n");
  for (int i=0; i<10; i++) printf(" "); printf("a_0");
  for (int i=0; i< 6; i++) printf(" "); printf("E_c");
  for (int i=0; i< 8; i++) printf(" "); printf("a_0");
  for (int i=0; i< 6; i++) printf(" "); printf("E_c"); printf("\n");
  for (int i=0; i<44; i++) printf("-"); printf("\n");

  for (int ii = 0; ii< ntype; ii++){
    printf("%6s ", elements[ii]);
    double amin = 1.5, amax = 4.5;
    double a0 = (amin + amax) * 0.5;
    const int nit_max = 1000;
    double Et;
    int n = 0;
    while (n < nit_max){
      double deda;
      lateng(a0, &Et, &deda, ii, 1);
      double da = deda * 0.05;
      double da2= da * da;
      if (da2 > 0.04) da /= sqrt(da2) * 5.;
      if (da2 < 1.e-12) break; 

      a0 -= da;
      if (a0 < amin) a0 = amin;
      if (a0 > amax) a0 = amax;
      //printf("n=%d, a0 = %g, da= %g\n", n, a0, da);
      n++;
    }
    //printf("%7.4f  %8.4f  ", a0, Et);
    if (n < nit_max) printf("%8.5f %8.4f  ", a0, Et);
    else printf("   ----     ----   ");

    amin = 2., amax = 6.;
    n = 0;
    a0 = (amin + amax) * 0.5;
    while (n < nit_max){
      double deda;
      lateng(a0, &Et, &deda, ii, 2);
      double da = deda * 0.05;
      double da2 = da * da;
      if (da2 > 0.04) da /= sqrt(da2) * 5.;
      if (da2 < 1.e-12) break;

      a0 -= da;
      if (a0 < amin) a0 = amin;
      if (a0 > amax) a0 = amax;
      //printf("n=%d, a0 = %g, da= %g\n", n, a0, da);
      n++;
    }
    if (n < nit_max) printf(" %8.5f %8.4f", a0, Et);
    else printf("   ----     ----");
    printf("\n");
  }
  for (int i=0; i<44; i++) printf("="); printf("\n");
}

void EAM::lateng(double a0, double *Epot, double *dEda, int ip, int lp)
{
  double Rshell[6], nShell[6];
  if (lp == 1){
    Rshell[0] = sqrt(3.)*0.5;
    Rshell[1] = 1.;
    Rshell[2] = sqrt(2.);
    Rshell[3] = sqrt(11.)*0.5;
    Rshell[4] = sqrt(3.);
    Rshell[5] = 2.;
  
    nShell[0] = 8.;
    nShell[1] = 6.;
    nShell[2] = 12;
    nShell[3] = 24.;
    nShell[4] = 8.;
    nShell[5] = 6.;
  } else {
    Rshell[0] = sqrt(2.)*0.5;
    Rshell[1] = 1.;
    Rshell[2] = sqrt(6.)*0.5;
    Rshell[3] = sqrt(2.);
    Rshell[4] = sqrt(10.)*0.5;
    Rshell[5] = sqrt(3.);

    nShell[0] = 12.;
    nShell[1] = 6.;
    nShell[2] = 24;
    nShell[3] = 12.;
    nShell[4] = 24.;
    nShell[5] = 8.;
  }

  *Epot = 0.;
  double rho = 0., rhop = 0., phip = 0.;

  for (int i=0; i<6; i++){
    double r = a0 * Rshell[i];
    if (r > rcut) break;
    rho  += Rho(r,  ip, ip) * nShell[i];
    rhop += Rhop(r, ip, ip) * Rshell[i] * nShell[i];
    phip += Phip(r, ip, ip) * Rshell[i] * nShell[i];
    *Epot += Phi(r, ip, ip) * nShell[i];
  }
  *Epot *= 0.5;
  *Epot += F(rho, ip);
  *dEda  = Fp(rho, ip) * rhop + 0.5 * phip;

return;
}
