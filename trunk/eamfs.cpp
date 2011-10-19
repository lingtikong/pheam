#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "eamfs.h"

using namespace std;

#define MAXLINE 256

/*------------------------------------------------------------------------------
 * Constructor is used to
 *----------------------------------------------------------------------------*/
EAMFS::EAMFS()
{
  char line[MAXLINE];
  // ask for the potential file
  do printf("\nPlease input the Setfl format EAM/FS potential file name: ");
  while (strlen(gets(line)) < 1);

  int n = strlen(line)+1;
  fname = new char[n];
  strcpy(fname, line);

  FILE *fp = fopen(fname,"r");
  if (fp == NULL) {
    printf("\nError: cannot open EAM potential file %s\n",fname);
    exit(1);
  }

  // read header, extract element names from nelements line
  fgets(line,MAXLINE,fp);
  fgets(line,MAXLINE,fp);
  fgets(line,MAXLINE,fp);
  fgets(line,MAXLINE,fp);

  sscanf(line,"%d", &ntype);

  int nwords = count_words(line);
  if (nwords != ntype+1){
    fclose(fp);
    printf("\nIncorrect element names in EAM potential file!\n");
    exit(2);
  }
  
  char **words = new char*[ntype+1];
  nwords = 0;
  char *first = strtok(line," \t\n\r\f");
  while (words[nwords++] = strtok(NULL," \t\n\r\f")) continue;

  elements = memory->create(elements,ntype, 10, "EAMFS_EAMFS:elements");
  for (int i=0; i<ntype; i++){
    strcpy(elements[i],words[i]);
  }
  delete [] words;

  fgets(line,MAXLINE,fp);
  Nrho = atoi(strtok(line," \t\n\r\f"));
  drho = atof(strtok(NULL," \t\n\r\f"));
  Nr   = atoi(strtok(NULL," \t\n\r\f"));
  dr   = atof(strtok(NULL," \t\n\r\f"));
  rcut = atof(strtok(NULL," \t\n\r\f"));

  mass = memory->create(mass, ntype, "EAMFS_EAMFS:mass");
  embed = memory->create(embed,ntype, Nrho, "EAMFS_EAMFS:embed");
  den   = memory->create(den, ntype, ntype, Nr, "EAMFS_EAMFS:den");
  pair  = memory->create(pair,ntype, ntype, Nr, "EAMFS_EAMFS:pair");

  for (int i = 0; i < ntype; i++) {
    fgets(line,MAXLINE,fp);
    strtok(line," \t\n\r\f");
    mass[i] = atof(strtok(NULL," \t\n\r\f"));

    grab(fp, Nrho, embed[i]);

    for (int j=0; j < ntype; j++) grab(fp, Nr, den[i][j]);
  }

  for (int i=0; i<ntype; i++){
    for (int j=0; j<=i; j++){
      grab(fp, Nr, pair[i][j]);
    }
  }
  fclose(fp);
  
  for (int i=0; i<ntype-1; i++){
    for (int j=i+1; j<ntype; j++){
      for (int k=0; k<Nr; k++) pair[i][j][k] = pair[j][i][k];
    }
  }
  
  // initialized other variables
  rdr   = 1./dr;
  rdrho = 1./drho;
  rcutsq = rcut*rcut;

  RHO0 = memory->create(RHO0, Nrho, "Setfl_Setfl:RHO0");
  R0   = memory->create(R0, Nr, "Setfl_Setfl:R0");
  for (int i=0; i<Nrho; i++) RHO0[i] = double(i)*drho;
  for (int i=0; i<Nr; i++) R0[i] = double(i)*dr;

  prepare_spline();

  initialized = 1;
}

/*------------------------------------------------------------------------------
 * Deconstructor to free memory
 *----------------------------------------------------------------------------*/
EAMFS::~EAMFS()
{
  for (int ip=0; ip<ntype; ip++){
    gsl_spline_free (F_spline[ip]);
    gsl_interp_accel_free (FAcc[ip]);
  }
  for (int jp=0; jp<ntype*ntype; jp++){
    gsl_spline_free (Den_spline[jp]);
    gsl_interp_accel_free (DenAcc[jp]);

    gsl_spline_free (Phi_spline[jp]);
    gsl_interp_accel_free (PhiAcc[jp]);
  }
  memory->sfree(FAcc);
  memory->sfree(DenAcc);
  memory->sfree(PhiAcc);
  memory->sfree(F_spline);
  memory->sfree(Den_spline);
  memory->sfree(Phi_spline);

  if (fname) delete []fname;
  memory->destroy(embed);
  memory->destroy(den);
  memory->destroy(pair);

  memory->destroy(R0);
  memory->destroy(RHO0);
}

/*------------------------------------------------------------------------------
 * F(rho)
 *----------------------------------------------------------------------------*/
double EAMFS::F(double rho, int ip)
{
  return gsl_spline_eval(F_spline[ip], rho, FAcc[ip]);
}

/*------------------------------------------------------------------------------
 * F'(rho)
 *----------------------------------------------------------------------------*/
double EAMFS::Fp(double rho, int ip)
{
  return gsl_spline_eval_deriv(F_spline[ip], rho, FAcc[ip]);
}

/*------------------------------------------------------------------------------
 * F"(rho)
 *----------------------------------------------------------------------------*/
double EAMFS::Fpp(double rho, int ip)
{
  return gsl_spline_eval_deriv2(F_spline[ip], rho, FAcc[ip]);
}

/*------------------------------------------------------------------------------
 * f(r)
 *----------------------------------------------------------------------------*/
double EAMFS::Rho(double r, int jp, int ip)
{
  int idx = jp * ntype + ip;
  return gsl_spline_eval(Den_spline[idx], r, DenAcc[idx]);
}

/*------------------------------------------------------------------------------
 * f'(r)
 *----------------------------------------------------------------------------*/
double EAMFS::Rhop(double r, int jp, int ip)
{
  int idx = jp * ntype + ip;
  return gsl_spline_eval_deriv(Den_spline[idx], r, DenAcc[idx]);
}

/*------------------------------------------------------------------------------
 * f"(r)
 *----------------------------------------------------------------------------*/
double EAMFS::Rhopp(double r, int jp, int ip)
{
  int idx = jp * ntype + ip;
  return gsl_spline_eval_deriv2(Den_spline[idx], r, DenAcc[idx]);
}

/*------------------------------------------------------------------------------
 * phi(r)
 *----------------------------------------------------------------------------*/
double EAMFS::Phi(double r, int ip, int jp)
{
  int idx = ip*ntype+jp;
  double pr = gsl_spline_eval(Phi_spline[idx], r, PhiAcc[idx]);
  return pr/r;
}

/*------------------------------------------------------------------------------
 * phi'(r)
 *----------------------------------------------------------------------------*/
double EAMFS::Phip(double r, int ip, int jp)
{
  int idx = ip*ntype+jp;
  double ppr = gsl_spline_eval_deriv(Phi_spline[idx], r, PhiAcc[idx]);
  double pr  = Phi(r, ip, jp);
  return (ppr - pr)/r;
}

/*------------------------------------------------------------------------------
 * phi"(r)
 *----------------------------------------------------------------------------*/
double EAMFS::Phipp(double r, int ip, int jp)
{
  int idx = ip*ntype+jp;
  double ppr = gsl_spline_eval_deriv2(Phi_spline[idx], r, PhiAcc[idx]);
  double pr  = Phip(r, ip, jp);
  return (ppr - pr - pr)/r;
}

/*------------------------------------------------------------------------------
 * Prepare for the splines
 *----------------------------------------------------------------------------*/
void EAMFS::prepare_spline()
{
  FAcc = (gsl_interp_accel **) memory->smalloc(ntype*sizeof(gsl_interp_accel *), "FAcc");
  F_spline = (gsl_spline **) memory->smalloc(ntype*sizeof(gsl_spline *), "F_spline");
  for (int ip=0; ip<ntype; ip++){
    FAcc[ip] = gsl_interp_accel_alloc ();
    F_spline[ip] = gsl_spline_alloc (gsl_interp_cspline, Nrho);

    gsl_spline_init(F_spline[ip], RHO0, embed[ip], Nrho);
  }

  DenAcc = (gsl_interp_accel **) memory->smalloc(ntype*ntype*sizeof(gsl_interp_accel *), "DenAcc");
  Den_spline = (gsl_spline **) memory->smalloc(ntype*ntype*sizeof(gsl_spline *), "Den_spline");
  int idx = 0;
  for (int ip=0; ip<ntype; ip++)
  for (int jp=0; jp<ntype; jp++){
    DenAcc[idx] = gsl_interp_accel_alloc ();
    Den_spline[idx] = gsl_spline_alloc (gsl_interp_cspline, Nr);

    gsl_spline_init(Den_spline[idx], R0, den[ip][jp], Nr);
    idx++;
  }

  PhiAcc = (gsl_interp_accel **) memory->smalloc(ntype*ntype*sizeof(gsl_interp_accel *), "PhiAcc");
  Phi_spline = (gsl_spline **) memory->smalloc(ntype*ntype*sizeof(gsl_spline *), "Phi_spline");
  idx=0;
  for (int ip=0; ip<ntype; ip++)
  for (int jp=0; jp<ntype; jp++){
    PhiAcc[idx] = gsl_interp_accel_alloc ();
    Phi_spline[idx] = gsl_spline_alloc (gsl_interp_cspline, Nr);

    gsl_spline_init(Phi_spline[idx], R0, pair[ip][jp], Nr);
    idx++;
  }

return;
}
/*----------------------------------------------------------------------------*/
