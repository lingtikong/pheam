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
    printf("  6. Local thermal properties by eigenvectors;\n");
    printf("  7. Phonon LDOS by Green's function;\n");
    printf("  0. Exit;\n");
    printf("Your choice [%d]: ", job);
    fgets(str,MAXLINE,stdin);
    char *ptr = strtok(str, " \t\n\r\f");
    if (ptr) job = atoi(ptr);
    for (int i=0; i<60; i++) printf("=");printf("\n");
  
    if (job < 1 || job >7) break;
    if (kpoints) delete kpoints;
    if (eig) dynmat->memory->destroy(eig);
    kpoints = NULL; eig = NULL;

    // Now to dispatch the job
    if (job <= 6){
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
          pldos(0);
          break;
        case 6:
          pldos(1);
          break;
        }
    } else if (job == 7){
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
void PHONON::pldos(int flag)
{
  char str[MAXLINE], id4ldos[MAXLINE], *ptr;
  //printf("\nThe # of atoms in cell is: %d, please input the atom IDs to compute\n", dynmat->natom);
  //printf("local PDOS, ID begins with 1. ");
  printf("\nThe # of atoms in cell is: %d, please input the atom IDs to compute\n", dynmat->natom);
  printf("local PDOS, IDs begin with 1. The input IDs follows an operator, which is\n");
  printf("one of =, >, <, >=, <=, <>, ><. For \"=\", one or more IDs can be appended;\n");
  printf("while for \">\", \">=\", \"<\" and \"<=\", only one is needed. For both \"<>\" and\n");
  printf("\"><\" two numbers are needed. Now please input you choice: ");
  int nword = dynmat->eam->count_words(fgets(str,MAXLINE,stdin));
  if ( nword < 2 ) return;
  nlocal = 0;
  ptr = strtok(str, " \t\n\r\f");
  if ( strcmp(ptr, "=") == 0 ){
     int nmax = nword-1;
     if (nmax > 0){
       locals = dynmat->memory->create(locals, nmax, "pldos:locals");
       strcpy(id4ldos, ptr);

       ptr = strtok(NULL, " \t\n\r\f");
       while (ptr != NULL){
         int id = atoi(ptr)-1;
         if (id >= 0 && id < dynmat->natom){
           int flag = 1;
           for (int i=0; i<nlocal; i++){
             if (id == locals[i]){ flag = 0; break; }
           }
           if (flag){
             locals[nlocal++] = id;
             strcat(id4ldos, " "); strcat(id4ldos, ptr);
           }
         }
         ptr = strtok(NULL, " \t\n\r\f");
       }
     }

  } else if (strcmp(ptr, ">") == 0){
     int nlow = atoi(strtok(NULL, " \t\n\r\f"));
     nlocal = dynmat->natom - nlow;
     if (nlocal < 1 || nlow < 0) return;

     locals = dynmat->memory->create(locals, nlocal, "pldos:locals");
     for (int i=0; i<nlocal; i++) locals[i] = nlow + i;
     sprintf(id4ldos, "> %d", nlow);

  } else if (strcmp(ptr, ">=") == 0){ 
     int nlow = atoi(strtok(NULL, " \t\n\r\f"));
     nlocal = dynmat->natom - nlow +1;
     if (nlocal < 1 || nlow < 1) return;

     locals = dynmat->memory->create(locals, nlocal, "pldos:locals");
     for (int i=0; i<nlocal; i++) locals[i] = nlow + i -1;
     sprintf(id4ldos, ">= %d", nlow);

  } else if (strcmp(ptr, "<") == 0){
     nlocal = atoi(strtok(NULL, " \t\n\r\f")) - 1;
     if (nlocal > dynmat->natom || nlocal < 1) return;

     locals = dynmat->memory->create(locals, nlocal, "pldos:locals");
     for (int i=0; i<nlocal; i++) locals[i] = i;

     sprintf(id4ldos, "< %d", nlocal+1);

  } else if (strcmp(ptr, "<=") == 0){
     nlocal = atoi(strtok(NULL, " \t\n\r\f"));
     if (nlocal > dynmat->natom || nlocal < 1) return;

     locals = dynmat->memory->create(locals, nlocal, "pldos:locals");
     for (int i=0; i<nlocal; i++) locals[i] = i;

     sprintf(id4ldos, "<= %d", nlocal);

  } else if (strcmp(ptr, "<>") == 0){
     int nlo = atoi(strtok(NULL, " \t\n\r\f")) -1;
     int nhi = atoi(strtok(NULL, " \t\n\r\f")) -1;
     if (nlo<0 || nhi >= dynmat->natom || nhi < nlo) return;

     nlocal = nhi - nlo + 1;
     locals = dynmat->memory->create(locals, nlocal, "pldos:locals");
     for (int i=0; i<nlocal; i++) locals[i] = i+nlo;

     sprintf(id4ldos, "[%d %d]", nlo+1, nhi+1);

  } else if (strcmp(ptr, "><") == 0){
     int nlo = atoi(strtok(NULL, " \t\n\r\f")) -1;
     int nhi = atoi(strtok(NULL, " \t\n\r\f")) -1;
     if (nlo<0 || nhi >= dynmat->natom || nhi < nlo) return;

     nlocal = nlo+1 + dynmat->natom - nhi;
     locals = dynmat->memory->create(locals, nlocal, "pldos:locals");
     for (int i=0; i<=nlo; i++) locals[i] = i;
     for (int i=nlo+1; i<nlocal; i++) locals[i] = nhi + i-(nlo+1);

     sprintf(id4ldos, "[1, %d] & [%d, %d]", nlo+1, nhi+1, dynmat->natom);

  } else {
     return;
  }

  if (nlocal < 1) return;
  printf("Local PDOS for %d atoms with id %s will be computed!\n", nlocal, id4ldos);

  const double tpi = 8.*atan(1.);
  wmin = 0.; wmax = 20.;
  printf("Please input the freqency (nv, THz) range to compute PDOS [%g %g]: ", wmin, wmax);
  if (dynmat->eam->count_words(fgets(str,MAXLINE,stdin)) >= 2){
    wmin = atof(strtok(str, " \t\n\r\f"));
    wmax = atof(strtok(NULL," \t\n\r\f"));
  }
  if (wmax < 0. || wmax < wmin) return;

  wmin *= tpi; wmax *= tpi;

  ndos = 201;
  printf("Please input your desired # of points in PDOS [%d]: ", ndos);
  fgets(str,MAXLINE,stdin);
  ptr = strtok(str," \t\n\r\f");
  if (ptr) ndos = atoi(ptr);

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
            double dr = dynmat->dm[idim][ipos+k].REALPART;
            double di = dynmat->dm[idim][ipos+k].IMAGPART;
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

  if (flag) compute_local_therm();

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
  if (dynmat->eam->count_words(fgets(str,MAXLINE,stdin)) >= 2){
    wmin = atof(strtok(str, " \t\n\r\f"));
    wmax = atof(strtok(NULL," \t\n\r\f"));
  }
  printf("You chose to get the DOS within [%lg %lg].\n", wmin, wmax);

  // initialize DOS variables
  ndos = 201;
  printf("Please input the number of points in the DOS [%d]: ", ndos);
  fgets(str,MAXLINE,stdin);
  char *ptr = strtok(str," \t\n\r\f");
  if (ptr) ndos = atoi(ptr);

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
  fgets(str,MAXLINE,stdin);
  char *ptr = strtok(str, " \t\n\r\f");
  if (ptr){
    int n = strlen(ptr)+1;
    fname = new char[n];
    strcpy(fname, ptr);
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

    fprintf(fp,"#Local phonon DOS for atom %d by eigenvalues\n", id);
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
  fgets(str,MAXLINE,stdin);
  char *ptr = strtok(str, " \t\n\r\f");
  if (ptr){
    int n = strlen(ptr) + 1;
    fname = new char[n];
    strcpy(fname, ptr);
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
  fgets(str,MAXLINE,stdin);
  char *ptr = strtok(str, " \t\n\r\f");
  if (ptr){
    int n = strlen(ptr)+1;
    fname = new char[n];
    strcpy(fname, ptr);
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
  fgets(str,MAXLINE,stdin); ptr = strtok(str, " \t\n\r\f");
  if (ptr) T = atof(ptr);

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
    printf("Done! Total cpu time used: %g second; wall time: %g seconds.\n", time->cpu_time(), time->wall_time());

    // ask for next temperature
    T = 0.;
    printf("Please input the desired temperature (K), non-positive to exit [0]: ");
    fgets(str,MAXLINE,stdin);
    char *ptr = strtok(str," \t\n\r\f"); if (ptr) T = atof(ptr);

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

/*------------------------------------------------------------------------------
 * Private method to compute the themal properties
 *----------------------------------------------------------------------------*/
void PHONON::compute_local_therm()
{
  char str[MAXLINE], *prefix;
  printf("\nPlease input the prefix for output local thermal files [lpth]: ");
  fgets(str,MAXLINE,stdin);
  char *ptr = strtok(str, " \t\n\r\f");
  if (ptr) {
    prefix = new char [strlen(ptr)+1];
    strcpy(prefix, ptr);
  } else {
    prefix = new char[5];
    strcpy(prefix, "lpth");
  }
  double Temp = 300.;
  printf("Please input the temperature to evaluate local thermal info [300]: ");
  fgets(str,MAXLINE,stdin);
  ptr = strtok(str, " \t\n\r\f");
  if (ptr) Temp = atof(ptr);

  Timer *time = new Timer();

  while ( Temp > 0. ){
    time->start();

    sprintf(str,"%s_%d.dat", prefix, int(Temp));
    FILE *fp = fopen(str, "w");
    fprintf(fp, "# Local thermal properties at %g K\n", Temp);
    fprintf(fp, "# 1    2  3  4  5  6  7  8  9  10 11 12 13 14 15 16 17\n");
    fprintf(fp, "# atom Ux Uy Uz Ut Sx Sy Sz St Fx Fy Fz Ft Cx Cy Cz Ct\n");
    fprintf(fp, "#      eV          kB          eV          kB\n");

    double h_o_KbT   = 6.6260755e1 / ( 1.380658e0 * Temp );
    double KbT_in_eV = 1.380658e-4 * Temp / 1.60217733e0;

    double Uvib[3], Svib[3], Fvib[3], Cvib[3];

    for (int ilocal = 0; ilocal<nlocal; ilocal++){
      int id = locals[ilocal];

      Uvib[0] = Uvib[1] = Uvib[2] = 0.;
      Svib[0] = Svib[1] = Svib[2] = 0.;
      Fvib[0] = Fvib[1] = Fvib[2] = 0.;
      Cvib[0] = Cvib[1] = Cvib[2] = 0.;

      double f = wmin-dw;
      for (int i=0; i<ndos; i++){
        f += dw;
        if (f <= 0.) continue;

        double x = h_o_KbT * f;
        double expx   = exp(x);
        double expxm1 = 1./(expx - 1.);

        double fac;
        if (i == 0) fac = 1.;
        else fac = double( (i%2+1)*2 );

        for (int idim=0; idim<3; idim++){
          Uvib[idim] += fac * ldos[i][ilocal][idim] * (0.5 + expxm1) * x;
          Svib[idim] += fac * ldos[i][ilocal][idim] * ( x * expxm1 - log(1.-1./expx));
          Fvib[idim] += fac * ldos[i][ilocal][idim] * log( 2.*sinh(0.5*x));
          Cvib[idim] += fac * ldos[i][ilocal][idim] * x * x * expx * expxm1 * expxm1;
        }
      }

      double dw3rd = dw / 3.;
      for (int idim=0; idim<3; idim++){
        Uvib[idim] *= dw3rd * KbT_in_eV;
        Svib[idim] *= dw3rd;
        Fvib[idim] *= dw3rd * KbT_in_eV;
        Cvib[idim] *= dw3rd;
      }

      fprintf(fp,"%d", id);
      for (int idim=0; idim<3; idim++) fprintf(fp," %lg", Uvib[idim]); fprintf(fp," %lg", Uvib[0]+Uvib[1]+Uvib[2]);
      for (int idim=0; idim<3; idim++) fprintf(fp," %lg", Svib[idim]); fprintf(fp," %lg", Svib[0]+Svib[1]+Svib[2]);
      for (int idim=0; idim<3; idim++) fprintf(fp," %lg", Fvib[idim]); fprintf(fp," %lg", Fvib[0]+Fvib[1]+Fvib[2]);
      for (int idim=0; idim<3; idim++) fprintf(fp," %lg", Cvib[idim]); fprintf(fp," %lg", Cvib[0]+Cvib[1]+Cvib[2]);
      fprintf(fp,"\n");
    }

    fclose(fp);
    time->stop();
    printf("Done! Total time used: %g second; wall time: %g seconds.\n", time->cpu_time(), time->wall_time());

    printf("Please input the temperature to evaluate local thermal info [%g]: ", Temp);
    fgets(str,MAXLINE,stdin);
    ptr = strtok(str, " \t\n\r\f"); if (ptr) Temp = atof(ptr);

  }
return;
}
