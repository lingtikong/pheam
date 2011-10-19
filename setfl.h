#ifndef SETFL_H
#define SETFL_H

#include "eam.h"
#include "gsl/gsl_spline.h"

class SETFL : public EAM {
public:
  SETFL();
  ~SETFL();

  double F(double, int);
  double Rho(double, int, int);
  double Phi(double, int, int);

  double Fp(double, int);
  double Rhop(double, int, int);
  double Phip(double, int, int);

  double Fpp(double, int);
  double Rhopp(double, int, int);
  double Phipp(double, int, int);

private:
  char *fname;
  int Nrho, Nr, AtNum;
  double rdrho, rdr;

  double **embed, **den, ***pair;

  double *RHO0, *R0;

  gsl_interp_accel **FAcc;
  gsl_interp_accel **DenAcc;
  gsl_interp_accel **PhiAcc;
  gsl_spline **F_spline;
  gsl_spline **Den_spline;
  gsl_spline **Phi_spline;

  void prepare_spline();

};

#endif
