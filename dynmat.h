#ifndef DYNMAT_H
#define DYNMAT_H

#include "eams.h"
#include "cell.h"
#include "memory.h"
#include "green.h"
#include <complex>

#ifdef MKL
#include "mkl.h"
#define COMPLEX MKL_Complex16
#define DOUBLEREAL double
#define REALPART real
#define IMAGPART imag
#define INTEGER  MKL_INT

#else

extern "C"{
#include "f2c.h"
#include "clapack.h"
}
#define COMPLEX doublecomplex
#define DOUBLEREAL doublereal
#define REALPART r
#define IMAGPART i
#define INTEGER integer
#endif

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
  COMPLEX **dm;

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
