#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "funcfl.h"

using namespace std;

#define MAXLINE 1024

/*------------------------------------------------------------------------------
 * Constructor is used to
 *----------------------------------------------------------------------------*/
FUNCFL::FUNCFL()
{
  char line[MAXLINE], potfiles[MAXLINE];
  // ask for the potential file
  do printf("\nPlease input the Funcfl format EAM potential file name(s): ");
  while (count_words(gets(potfiles)) < 1);

  int nfile = count_words(potfiles);
  char **words = new char*[nfile+1];
  int nwords = 0;
  words[0] = strtok(potfiles," \t\n\r\f");
  while (words[++nwords] = strtok(NULL," \t\n\r\f")) continue;
  mass = memory->create(mass,nfile, "FUNCFL:mass");
  elements = memory->create(elements,nfile, 10, "Funcfl_Funcfl:elements");

  ntype = 0;
  for (int ii = 0; ii<nfile; ii++){
    int n = strlen(words[ii])+1;
    char * fname = new char[n];
    strcpy(fname, words[ii]);
    FILE *fp = fopen(fname,"r");

    if (fp == NULL) {
      printf("\nError: Cannot open EAM potential file %s.\n", fname);
      delete []fname;
      continue;
    }

    // read header, extract element names from nelements line
    fgets(line,MAXLINE,fp);
    fgets(line,MAXLINE,fp);
    strcpy(elements[ntype], strtok(line," \t\n\r\f"));
    mass[ntype] = atof(strtok(NULL," \t\n\r\f"));
  
    int Nrhop, Nrp;
    double drhop, drp, rcutp;
    fgets(line,MAXLINE,fp);
    Nrhop = atoi(strtok(line," \t\n\r\f"));
    drhop = atof(strtok(NULL," \t\n\r\f"));
    Nrp   = atoi(strtok(NULL," \t\n\r\f"));
    drp   = atof(strtok(NULL," \t\n\r\f"));
    rcutp = atof(strtok(NULL," \t\n\r\f"));

    if (ntype ==0){
      Nrho = Nrhop; drho = drhop; Nr = Nrp; dr = drp; rcut = rcutp;
      embed = memory->create(embed,nfile, Nrho, "Funcfl_Funcfl:embed");
      den   = memory->create(den,nfile, Nr, "Funcfl_Funcfl:den");
      pair  = memory->create(pair,nfile, Nr, "Funcfl_Funcfl:pair");
    } else {
      if (Nrho != Nrhop || Nr != Nrp || abs(drhop-drho)>1.e-9 || abs(drp-dr) > 1.e-9 || abs(rcutp-rcut)>1.e-5){
        printf("\nError: Header info read from %s is different from previous ones! Current potential ignored.\n", fname);
        delete []fname;
        continue;
      }
    }

    grab(fp,Nrho,embed[ntype]);
    grab(fp,Nr,   pair[ntype]);
    grab(fp,Nr,    den[ntype]);
    fclose(fp);
    
    delete []fname;

    ntype++;
  }
  delete [] words;

  RHO0 = memory->create(RHO0, Nrho, "Setfl_Setfl:RHO0");
  R0   = memory->create(R0, Nr, "Setfl_Setfl:R0");

  for (int i=0; i<Nrho; i++) RHO0[i] = double(i)*drho;
  for (int i=0; i<Nr; i++) R0[i] = double(i)*dr;

  // initialized other variables
  rdr   = 1./dr;
  rdrho = 1./drho;
  rcutsq = rcut*rcut;

  prepare_spline();

  initialized = 1;
}

/*------------------------------------------------------------------------------
 * Deconstructor to free memory
 *----------------------------------------------------------------------------*/
FUNCFL::~FUNCFL()
{
  for (int ip=0; ip<ntype; ip++){
    gsl_spline_free (F_spline[ip]);
    gsl_spline_free (Den_spline[ip]);
    gsl_interp_accel_free (FAcc[ip]);
    gsl_interp_accel_free (DenAcc[ip]);

    gsl_spline_free (Phi_spline[ip]);
    gsl_interp_accel_free (PhiAcc[ip]);
  }
  memory->sfree(FAcc);
  memory->sfree(DenAcc);
  memory->sfree(PhiAcc);
  memory->sfree(F_spline);
  memory->sfree(Den_spline);
  memory->sfree(Phi_spline);

  memory->destroy(embed);
  memory->destroy(den);
  memory->destroy(pair);

  memory->destroy(R0);
  memory->destroy(RHO0);

}

/*------------------------------------------------------------------------------
 * F(rho)
 *----------------------------------------------------------------------------*/
double FUNCFL::F(double rho, int ip)
{
  return gsl_spline_eval(F_spline[ip], rho, FAcc[ip]);
}

/*------------------------------------------------------------------------------
 * F'(rho)
 *----------------------------------------------------------------------------*/
double FUNCFL::Fp(double rho, int ip)
{
  return gsl_spline_eval_deriv(F_spline[ip], rho, FAcc[ip]);
}

/*------------------------------------------------------------------------------
 * F"(rho)
 *----------------------------------------------------------------------------*/
double FUNCFL::Fpp(double rho, int ip)
{
  return gsl_spline_eval_deriv2(F_spline[ip], rho, FAcc[ip]);
}

/*------------------------------------------------------------------------------
 * f(r)
 *----------------------------------------------------------------------------*/
double FUNCFL::Rho(double r, int jp, int ip)
{
  return gsl_spline_eval(Den_spline[jp], r, DenAcc[jp]);
}

/*------------------------------------------------------------------------------
 * f'(r)
 *----------------------------------------------------------------------------*/
double FUNCFL::Rhop(double r, int jp, int ip)
{
  return gsl_spline_eval_deriv(Den_spline[jp], r, DenAcc[jp]);
}

/*------------------------------------------------------------------------------
 * f"(r)
 *----------------------------------------------------------------------------*/
double FUNCFL::Rhopp(double r, int jp, int ip)
{
  return gsl_spline_eval_deriv2(Den_spline[jp], r, DenAcc[jp]);
}

/*------------------------------------------------------------------------------
 * phi(r)
 *----------------------------------------------------------------------------*/
double FUNCFL::Phi(double r, int ip, int jp)
{
  const double RydBohr2evA = 27.2*0.529;
  double pi, pj;
  pi = gsl_spline_eval(Phi_spline[ip], r, PhiAcc[ip]);
  if (ip == jp) pj = pi;
  else pj = gsl_spline_eval(Phi_spline[jp], r, PhiAcc[jp]);

  return RydBohr2evA * (pi*pj)/r;
}

/*------------------------------------------------------------------------------
 * phi'(r)
 *----------------------------------------------------------------------------*/
double FUNCFL::Phip(double r, int ip, int jp)
{
  const double RydBohr2evA = 27.2*0.529;
  double pi, pj, pip, pjp;
  pi  = gsl_spline_eval(Phi_spline[ip], r, PhiAcc[ip]);
  pip = gsl_spline_eval_deriv(Phi_spline[ip], r, PhiAcc[ip]);
  if (ip == jp){
    pj  = pi;
    pjp = pip;
  } else {
    pj  = gsl_spline_eval(Phi_spline[jp], r, PhiAcc[jp]);
    pjp = gsl_spline_eval_deriv(Phi_spline[jp], r, PhiAcc[jp]);
  }

  return RydBohr2evA * (pip*pj + pi*pjp - pi*pj/r)/r;
}

/*------------------------------------------------------------------------------
 * phi"(r)
 *----------------------------------------------------------------------------*/
double FUNCFL::Phipp(double r, int ip, int jp)
{
  const double RydBohr2evA = 27.2*0.529;
  double pi, pj, pip, pjp, pipp, pjpp;
  pi   = gsl_spline_eval(Phi_spline[ip], r, PhiAcc[ip]);
  pip  = gsl_spline_eval_deriv(Phi_spline[ip], r, PhiAcc[ip]);
  pipp = gsl_spline_eval_deriv2(Phi_spline[ip], r, PhiAcc[ip]);
  if (ip == jp){
    pj   = pi;
    pjp  = pip;
    pjpp = pipp;
  } else {
    pj   = gsl_spline_eval(Phi_spline[jp], r, PhiAcc[jp]);
    pjp  = gsl_spline_eval_deriv(Phi_spline[jp], r, PhiAcc[jp]);
    pjpp = gsl_spline_eval_deriv2(Phi_spline[jp], r, PhiAcc[jp]);
  }

  double r2 = r*r;

  return RydBohr2evA * (pipp*pj + 2.*pip*pjp - 2.*pip*pj/r + pi*pjpp - 2.*pi*pjp/r + 2.*pi*pj/r2) / r;
}

/*------------------------------------------------------------------------------
 * Prepare for the splines
 *----------------------------------------------------------------------------*/
void FUNCFL::prepare_spline()
{
  FAcc = (gsl_interp_accel **) memory->smalloc(ntype*sizeof(gsl_interp_accel *), "FAcc");
  F_spline = (gsl_spline **) memory->smalloc(ntype*sizeof(gsl_spline *), "F_spline");
  for (int ip=0; ip<ntype; ip++){
    FAcc[ip] = gsl_interp_accel_alloc ();
    F_spline[ip] = gsl_spline_alloc (gsl_interp_cspline, Nrho);

    gsl_spline_init(F_spline[ip], RHO0, embed[ip], Nrho);
  }

  DenAcc = (gsl_interp_accel **) memory->smalloc(ntype*sizeof(gsl_interp_accel *), "DenAcc");
  Den_spline = (gsl_spline **) memory->smalloc(ntype*sizeof(gsl_spline *), "Den_spline");
  for (int ip=0; ip<ntype; ip++){
    DenAcc[ip] = gsl_interp_accel_alloc ();
    Den_spline[ip] = gsl_spline_alloc (gsl_interp_cspline, Nr);

    gsl_spline_init(Den_spline[ip], R0, den[ip], Nr);
  }

  PhiAcc = (gsl_interp_accel **) memory->smalloc(ntype*sizeof(gsl_interp_accel *), "PhiAcc");
  Phi_spline = (gsl_spline **) memory->smalloc(ntype*sizeof(gsl_spline *), "Phi_spline");
  for (int ip=0; ip<ntype; ip++){
    PhiAcc[ip] = gsl_interp_accel_alloc ();
    Phi_spline[ip] = gsl_spline_alloc (gsl_interp_cspline, Nr);

    gsl_spline_init(Phi_spline[ip], R0, pair[ip], Nr);
  }

return;
}
/*----------------------------------------------------------------------------*/
