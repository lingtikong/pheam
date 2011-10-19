#ifndef PHONON_H
#define PHONON_H

#include "dynmat.h"
#include "kpoints.h"

class PHONON{
public:
  PHONON();
  ~PHONON();

  DYNMAT *dynmat;
  KPOINTS *kpoints;

private:
  void ComputeAll();

  void pdos();
  void pdisp();
  void pldos();
  void therm();

  void writedos();
  void writeldos();

  void normalize();

  int nq, ndim;
  double **q, *w, *qr;

  int ndos;
  double wmin, wmax, dw, rdw;
  int *locals, nlocal;
  double **eig;
  double *dos, ***ldos;
};
#endif
