#ifndef DYNMAT_H
#define DYNMAT_H

#include "eams.h"
#include "cell.h"
#include "memory.h"
#include "green.h"

extern "C"{
#include "f2c.h"
#include "clapack.h"
}

class DYNMAT {
public:
  DYNMAT();
  ~DYNMAT();

  CELL *atom;
  EAM  *eam;
  Green *green;
  Memory *memory;

  int natom, ndim;
  double Ec, Ep, Et;
  int *map;    // to map atomic type in cell to that in EAM

  void computeDM(double *);
  int  computeEigen(int); // int flag: 0, eigenvalues only; 1: both eigenvalues and eigenvectors
  void GreenLDOS();

  double *egval;
  doublecomplex **dm;

private:
  void selectEAM(void);
  void checkmap(void);
  void setup();

  int neimax, Lx, Ly, Lz;
  int *NumNei, **NeiList;
  double ***Bonds;
  double *den, *Fp, *Fpp;

};
#endif
