#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "timer.h"
#include "phonon.h"

#define MAX(a,b) ((a) > (b) ? (a) : (b))
#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAXLINE 256

/*------------------------------------------------------------------------------
 * Constructor is used as the driver of the phonon calculations based on EAM
 *----------------------------------------------------------------------------*/
PHONON::PHONON()
{
  dynmat = new DYNMAT;
  eig     = NULL;
  dos     = NULL;
  ldos    = NULL;
  locals  = NULL;
  kpoints = NULL;

  int job = 0;
  char str[MAXLINE];
  // main menu
  while (1){
    printf("\n");for (int i=0; i<60; i++) printf("=");printf("\n");
    printf("Please select the job to do:\n");
    printf("  1. Gamma-point phonon DOS calculation;\n");
    printf("  2. Full phonon DOS calculations;\n");
    printf("  3. Phonon dispersion calculations;\n");
    printf("  4. Thermal properties calculations;\n");
    printf("  5. Phonon LDOS by eigenvectors;\n");
    printf("  6. Phonon LDOS by Green's function;\n");
    printf("  0. Exit;\n");
    printf("Your choice [%d]: ", job);
    if (strlen(gets(str)) > 0) sscanf(str,"%d", &job);
    for (int i=0; i<60; i++) printf("=");printf("\n");
  
    if (job < 1 || job >6) break;
    if (kpoints) delete kpoints;
    if (eig) dynmat->memory->destroy(eig);
    kpoints = NULL; eig = NULL;

    // Now to dispatch the job
    if (job <= 5){
      kpoints = new KPOINTS(job, dynmat->atom);
      
      ndim = dynmat->ndim;
      nq = kpoints->nq;
      q  = kpoints->q;
      w  = kpoints->w;
      qr = kpoints->qr;
      eig  = dynmat->memory->create(eig,nq, ndim, "PHONON_PHONON:eig");
  
      if (job <= 4) ComputeAll();

      switch (job){
        case 1:
        case 2:
          pdos();
          break;
        case 3:
          pdisp();
          break;
        case 4:
          therm();
          break;
        case 5:
          pldos();
          break;
        }
    } else if (job == 6){
      dynmat->GreenLDOS();
    }
  }
return;
}

/*------------------------------------------------------------------------------
 * Deconstructor is used to free all memory
 *----------------------------------------------------------------------------*/
PHONON::~PHONON()
{
  dynmat->memory->destroy(eig);
  dynmat->memory->destroy(ldos);
  dynmat->memory->destroy(dos);
  dynmat->memory->destroy(locals);
  
  if (kpoints) delete kpoints;
  if (dynmat)  delete dynmat;

  q  = NULL;
  w  = NULL;
  qr = NULL;

}

/*------------------------------------------------------------------------------
 * Private method, to calculate the local phonon DOS based on the eigenvectors
 *----------------------------------------------------------------------------*/
void PHONON::pldos()
{
  char str[MAXLINE], *ptr;
  printf("\nThe # of atoms in cell is: %d, please input the atom IDs to compute\n", dynmat->natom);
  printf("local PDOS, ID begins with 1: ");
  int nmax = dynmat->eam->count_words(gets(str));
  if (nmax < 1) return;
  locals = dynmat->memory->create(locals, nmax, "pldos:locals");

  nlocal = 0;
  ptr = strtok(str," \t\n\r\f");
  while (ptr != NULL){
    int id = atoi(ptr)-1;
    if (id >= 0 && id < dynmat->natom) locals[nlocal++] = id;
  
    ptr = strtok(NULL," \t\n\r\f");
  }
  if (nlocal < 1) return;
  printf("Local PDOS for atom(s):");
  for (int i=0; i<nlocal; i++) printf(" %d", locals[i]+1);
  printf(" will be computed.\n");

  const double tpi = 8.*atan(1.);
  wmin = 0.; wmax = 10.;
  printf("Please input the freqency (nv, THz) range to compute PDOS [%g %g]: ", wmin, wmax);
  if (dynmat->eam->count_words(gets(str)) >= 2) sscanf(str,"%lg %lg", &wmin, &wmax);
  if (wmax < 0. || wmax < wmin) return;

  wmin *= tpi; wmax *= tpi;

  ndos = 101;
  printf("Please input your desired # of points in PDOS [%d]: ", ndos);
  if (strlen(gets(str)) > 0) ndos = atoi(strtok(str," \t\n\r\f"));
  if (ndos < 1) return;
  ndos += (ndos+1)%2;

  dw = (wmax-wmin)/double(ndos-1);
  rdw = 1./dw;

  ldos = dynmat->memory->create(ldos,ndos,nlocal,3,"phonon_pldos:ldos");
  dos  = dynmat->memory->create(dos,ndos,"phonon_pldos:dos");

  for (int i=0; i<ndos; i++){
    dos[i] = 0.;
    for (int j=0; j<nlocal; j++)
    for (int k=0; k<3; k++) ldos[i][j][k] = 0.;
  }

  int nprint;
  if (nq > 10) nprint = nq/10;
  else nprint = 1;

  Timer *time = new Timer();

  printf("\nNow to compute the phonons "); fflush(stdout);
  for (int iq=0; iq<nq; iq++){
    if ((iq+1)%nprint == 0) {printf("."); fflush(stdout);}

    dynmat-> computeDM(q[iq]);
    if (dynmat->computeEigen(1) != 0){
      printf("\nError while computing the eigen problem for the %d-th q: [%lg %lg %lg]\n",
      iq+1, q[iq][0], q[iq][1], q[iq][2]);
      exit(1);
    }
    for (int idim=0; idim<ndim; idim++){
      int hit = int((dynmat->egval[idim] - wmin)*rdw+0.5);
      if (hit >= 0 && hit <ndos){
        dos[hit] += w[iq];

        for (int ilocal=0; ilocal<nlocal; ilocal++){
          int ipos = locals[ilocal]*3; 
          for (int k=0; k<3; k++){
            double dr = dynmat->dm[idim][ipos+k].r, di = dynmat->dm[idim][ipos+k].i;
            double norm = dr * dr + di * di;
            ldos[hit][ilocal][k] += w[iq] * norm;
          }
        }
      }
    }
  }

  // convert omega to nu
  wmin /= tpi; wmax /=tpi; dw /= tpi; rdw *=tpi;

  normalize();
  time->stop(); time->print(); delete time;

  writedos();
  writeldos();

  dynmat->memory->destroy(ldos);   ldos = NULL;
  dynmat->memory->destroy(dos);    dos  = NULL;
  dynmat->memory->destroy(locals); locals = NULL;

return;
}

/*------------------------------------------------------------------------------
 * Private method, as the driver to compute all eigenvalues
 *----------------------------------------------------------------------------*/
void PHONON::ComputeAll()
{
  Timer *time = new Timer();

  int info, nprint;
  if (nq > 10) nprint = nq/10;
  else nprint = 1;
  printf("\nNow to compute the phonons "); fflush(stdout);

  for (int iq=0; iq<nq; iq++){
    if ((iq+1)%nprint == 0) {printf("."); fflush(stdout);}

    dynmat-> computeDM(q[iq]);
    if (dynmat->computeEigen(0) == 0){
      for (int i=0; i<ndim; i++) eig[iq][i] = dynmat->egval[i];
    } else {
      printf("\nError while computing the eigen problem for the %d-th q: [%lg %lg %lg]\n",
      iq+1, q[iq][0], q[iq][1], q[iq][2]);
      exit(1);
    }
  }

  time->stop(); time->print(); delete time;

return;
}

/*------------------------------------------------------------------------------
 * Private method to compute the phonon DOS
 *----------------------------------------------------------------------------*/
void PHONON::pdos()
{
  // determin the frequency range to compute DOS
  char str[MAXLINE];
  const double rtpi = 1./(8.*atan(1.));

  wmin = wmax = eig[0][0]*rtpi;
  for (int iq=0; iq<nq; iq++){
    for (int i=0; i<ndim; i++){eig[iq][i] *= rtpi; wmin=MIN(wmin, eig[iq][i]); wmax=MAX(wmax,eig[iq][i]);}
  }
  printf("\nPlease input the range to evaluate PDOS [%lg %lg]: ", wmin, wmax);
  if (dynmat->eam->count_words(gets(str)) >= 2) sscanf(str, "%lg %lg", &wmin, &wmax);
  printf("You chose to get the DOS within [%lg %lg].\n", wmin, wmax);

  // initialize DOS variables
  ndos = 201;
  printf("Please input the number of points in the DOS [%d]: ", ndos);
  if (strlen(gets(str)) > 0) ndos = atoi(strtok(str," \t\n\r\f"));
  if (ndos < 2) return;
  ndos += (ndos+1)%2;

  dos = dynmat->memory->create(dos,ndos,"pdos:dos");
  for (int i=0; i<ndos; i++) dos[i] = 0.;

  dw = (wmax-wmin)/double(ndos-1);
  rdw = 1./dw;

  Timer *time = new Timer();
  printf("\nNow to compute the phonon dos ...");
  // compute the DOS
  double offset = wmin-0.5*dw;
  for (int iq=0; iq<nq; iq++){
    for (int i=0; i<ndim; i++){
      int idos = int((eig[iq][i] - offset)*rdw);
      if (idos < ndos && idos >=0) dos[idos] += w[iq];
    }
  }

  // normalize the DOS
  normalize();

  printf("Done!\n");
  time->stop(); time->print(); delete time;

  // output the DOS
  writedos();

  dynmat->memory->destroy(dos); dos = NULL;

return;
}

/*------------------------------------------------------------------------------
 * Private method to write the pnonon DOS to file
 *----------------------------------------------------------------------------*/
void PHONON::writedos()
{
  char str[MAXLINE], *fname;
  printf("\nPlease input the file name to output the phonon DOS [pdos.dat]: ");
  if (strlen(gets(str)) > 1){
    int n = strlen(str)+1;
    fname = new char[n];
    strcpy(fname, str);
  } else {
    fname = new char[9];
    strcpy(fname,"pdos.dat");
  }
  FILE *fp = fopen(fname, "w");

  fprintf(fp,"#nu DOS\n");
  double wnow = wmin;
  for (int i=0; i<ndos; i++){
    fprintf(fp,"%lg %lg\n",wnow, dos[i]);
    wnow += dw;
  }
  fclose(fp);

  delete []fname;

return;
}

/*------------------------------------------------------------------------------
 * Private method to write the pnonon local DOS to file
 *----------------------------------------------------------------------------*/
void PHONON::writeldos()
{
  char str[MAXLINE];
  for (int ilocal=0; ilocal<nlocal; ilocal++){
    int id = locals[ilocal]+1;
    int n = sprintf(str,"pldos_%d.dat", id);
    char *fname = new char[n+1];
    strcpy(fname, str);
    FILE *fp = fopen(fname, "w");

    fprintf(fp,"#Local phonon DOS for atom %d\n", id);
    fprintf(fp,"#freq x  y  z  total\n");
    
    double wnow = wmin;
    for (int idos=0; idos<ndos; idos++){
      double ldx = ldos[idos][ilocal][0];
      double ldy = ldos[idos][ilocal][1];
      double ldz = ldos[idos][ilocal][2];
      double total = ldx + ldy + ldz;
      fprintf(fp,"%lg %lg %lg %lg %lg\n", wnow, ldx, ldy, ldz, total);
      wnow += dw;
    }
    fclose(fp);
    delete []fname;
  }

return;
}

/*------------------------------------------------------------------------------
 * Private method to compute the phonon dispersion curve
 *----------------------------------------------------------------------------*/
void PHONON::pdisp()
{
  // ask for the file name
  char str[MAXLINE], *fname;
  printf("\nPlease input the file name to output the phonon dispersion [pdisp.dat]: ");
  if (strlen(gets(str)) > 0){
    int n = strlen(str) + 1;
    fname = new char[n];
    strcpy(fname, str);
  } else {
    fname = new char[10];
    strcpy(fname,"pdisp.dat");
  }
  FILE *fp = fopen(fname, "w");
  
  const double rtpi = 1./(8.*atan(1.));
  fprintf(fp,"#index qr qx qy qz %d-frequencies\n", ndim);

  for (int iq=0; iq<nq; iq++){
    fprintf(fp,"%d %lg %lg %lg %lg", iq+1, qr[iq], q[iq][0], q[iq][1], q[iq][2]);
    for (int i=0; i<ndim; i++) fprintf(fp," %lg", eig[iq][i]*rtpi); fprintf(fp,"\n");
  }
  fclose(fp);

  delete []fname;

return;
}

/*------------------------------------------------------------------------------
 * Private method to compute the themal properties
 *----------------------------------------------------------------------------*/
void PHONON::therm()
{
  // ask for file name to output result
  char str[MAXLINE], *fname;
  printf("\nPlease input the file name to output thermal properties [pthermo.dat]: ");
  if (strlen(gets(str)) > 0){
    int n = strlen(str)+1;
    fname = new char[n];
    strcpy(fname, str);
  } else {
    fname = new char[12];
    strcpy(fname,"pthermo.dat");
  }
  FILE *fp = fopen(fname, "w");
  // header line 
  fprintf(fp,"#Temp  Uvib Svib Fvib ZPE Cvib\n");
  fprintf(fp,"# K    eV   Kb   eV   eV  Kb\n");

  // constants
  const double hbar = 1.05457266e-34, Kb = 1.380658e-23, eV = 1.60217733e-19;

  // first temperature
  double T = 0.;
  printf("Please input the desired temperature (K), non-positive to exit [0]: ");
  if (strlen(gets(str)) > 0) T = atof(strtok(str," \t\n\r\f"));

  Timer *time = new Timer();

  while (T > 0.){
    time->start();
    printf("Now to compute for temperature %g ... ", T);
    // constants under the same temperature; assuming angular frequency in THz
    double hb_o_KbT = hbar/(Kb*T)*1.e12, KbT_in_eV = Kb*T/eV;
    
    double Uvib = 0., Svib = 0., Fvib = 0., Cvib = 0., ZPE = 0.;
    for (int iq=0; iq<nq; iq++){
      double Utmp = 0., Stmp = 0., Ftmp = 0., Ztmp = 0., Ctmp = 0.;
      for (int i=0; i<ndim; i++){
        if (eig[iq][i] <= 0.) continue;
        double x = eig[iq][i] * hb_o_KbT;
        double expterm = 1./(exp(x)-1.);
        Stmp += x*expterm - log(1.-exp(-x));
        Utmp += (0.5+expterm)*x;
        Ftmp += log(2.*sinh(0.5*x));
        Ctmp += x*x*exp(x)*expterm*expterm;
        Ztmp += 0.5*hbar*eig[iq][i];
      }
      Svib += w[iq]*Stmp;
      Uvib += w[iq]*Utmp;
      Fvib += w[iq]*Ftmp;
      Cvib += w[iq]*Ctmp;
      ZPE  += w[iq]*Ztmp;
    }
    Uvib *= KbT_in_eV;
    Fvib *= KbT_in_eV;
    ZPE  /= eV*1.e-12;
    // output result under current temperature
    fprintf(fp,"%lg %lg %lg %lg %lg %lg\n", T, Uvib, Svib, Fvib, ZPE, Cvib);

    time->stop();
    printf("Done! Total time used: %g second.\n", time->elapse());
    // ask for next temperature
    printf("Please input the desired temperature (K), non-positive to exit [0]: ");
    if (strlen(gets(str)) > 0) T = atof(strtok(str," \t\n\r\f"));
    else T = 0.;
  }
  fclose(fp);
  delete time;

return;
}

/*------------------------------------------------------------------------------
 * Private method, to normalize the calcated phonon dos and/or local dos
 *----------------------------------------------------------------------------*/
void PHONON::normalize()
{
  if (dos){ // to normalize dos
    double odd = 0., even = 0.;
    for (int i=1; i<ndos-1; i+=2) odd  += dos[i];
    for (int i=2; i<ndos-1; i+=2) even += dos[i];

    double sum = dos[0] + dos[ndos-1];
    sum += 4.*odd+2.*even;
    sum = 3.*rdw/sum;
    for (int i=0; i<ndos; i++) dos[i] *= sum;
  }

  if (ldos){ // to normalize ldos
    for (int j=0; j<nlocal; j++)
    for (int k=0; k<3; k++){
      double odd = 0., even = 0.;
      for (int i=1; i<ndos-1; i+=2) odd  += ldos[i][j][k];
      for (int i=2; i<ndos-1; i+=2) even += ldos[i][j][k];

      double sum = ldos[0][j][k] + ldos[ndos-1][j][k];
      sum += 4.*odd+2.*even;
      sum = 3.*rdw/sum;
      for (int i=0; i<ndos; i++) ldos[i][j][k] *= sum;
    }
  }
return;
}
