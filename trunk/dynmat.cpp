#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "timer.h"
#include "green.h"

#include "dynmat.h"
#ifdef OMP
#include "omp.h"
#endif

using namespace std;

#define MAX(a,b) ((a) > (b) ? (a) : (b))
#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAXLINE 256

/*------------------------------------------------------------------------------
 * Constructor is used to allocate memory and initialize necessary variables
 *----------------------------------------------------------------------------*/
DYNMAT::DYNMAT()
{
  // create objects
  memory = new Memory;
  atom   = new CELL;
  selectEAM();

  // get the mapping between atomic types from "atom" and those from "EAM"
  memory->create(map,atom->ntype,"DYNMAT:map");
  checkmap();

  // allocate variables
  natom = atom->natom;
  ndim = natom*3;
  neimax = 40;

  memory->create(den, natom, "DYNMAT:den");
  memory->create(Fp , natom, "DYNMAT:Fp");
  memory->create(Fpp, natom, "DYNMAT:Fpp");
  memory->create(NumNei, natom, "DYNMAT:NumNei");

  memory->create(NeiList,neimax, natom,"DYNMAT_DYNMAT:NeiList");
  memory->create(Bonds,neimax,natom,4,"DYNMAT_DYNMAT:Bonds");

  // ask for the extension of lattice in three dimensions
  // assuming orthogonal lattice
  Lx = int(eam->rcut/atom->axis[0][0]+1.0); if (atom->pbc[0] == 0) Lx=0;
  Ly = int(eam->rcut/atom->axis[1][1]+1.0); if (atom->pbc[1] == 0) Ly=0;
  Lz = int(eam->rcut/atom->axis[2][2]+1.0); if (atom->pbc[2] == 0) Lz=0;
  char str[MAXLINE];
  printf("\nPlease input the # of extension of cell in three dimensions [%d %d %d]: ", Lx,Ly,Lz);
  if (eam->count_words(fgets(str,MAXLINE,stdin)) >= 3){
    Lx = atof(strtok(str, " \t\n\r\f"));
    Ly = atof(strtok(NULL," \t\n\r\f"));
    Lz = atof(strtok(NULL," \t\n\r\f"));
  }
  if (atom->pbc[0] == 0) Lx=0;
  if (atom->pbc[1] == 0) Ly=0;
  if (atom->pbc[2] == 0) Lz=0;
  printf("The # of extension in three dimensions accepted: %d x %d x %d.", Lx,Ly,Lz);

  // now to get the energy, electronic density and embedded terms of the system
  setup();

  // allocate memory for dynamical matrix and its eigen values/vectors
  memory->create(egval,ndim,"DYNMAT_DYNMAT:egval");
  memory->create(dm,ndim, ndim,"DYNMAT_DYNMAT:dm");

return;
}

/*------------------------------------------------------------------------------
 * Deconstructor is used to free memory
 *----------------------------------------------------------------------------*/
DYNMAT::~DYNMAT()
{
  memory->destroy(map);
  memory->destroy(den);
  memory->destroy(Fp);
  memory->destroy(Fpp);
  memory->destroy(NumNei);
  memory->destroy(egval);
  memory->destroy(NeiList);
  memory->destroy(Bonds);
  memory->destroy(dm);

  delete eam;
  delete atom;
  delete memory;
}

/*------------------------------------------------------------------------------
 * Private method to get the mapping between atomic types from "atom" and
 * those from "eam"
 *----------------------------------------------------------------------------*/
void DYNMAT::checkmap()
{
  char str[MAXLINE];
  for (int i=0; i<atom->ntype; i++){
    int ip = eam->index(atom->elements[i]);
    if (ip < 0){
      printf("\nWarning: Cannot find element %s in EAM; available elements in EAM are:", atom->elements[i]);
      for (int j=0; j<eam->ntype; j++) printf(" %s", eam->elements[j]); printf("\n");
      while (ip<0){
        char ename[10];
        printf("Please input the EAM name for %s [%s]: ", atom->elements[i], eam->elements[0]);
        if (eam->count_words(fgets(str,MAXLINE,stdin)) < 1)
          strcpy(ename, eam->elements[0]);
        else strcpy(ename, strtok(str, " \n\t\r\f"));

        ip = eam->index(ename);
      }
    }
    map[i] = ip;
  }
return;
}

/*------------------------------------------------------------------------------
 * Private method to setup the calculation
 *----------------------------------------------------------------------------*/
void DYNMAT::setup()
{
  Ec = Ep = 0.;
  // get the electronic density and neighbors
#ifdef OMP
  int npmax = omp_get_max_threads();
  omp_set_num_threads(npmax);
#endif
  #pragma omp parallel for default(shared) schedule(guided)
  for (int i=0; i< natom; i++) den[i] = 0.;

  #pragma omp parallel for default(shared) schedule(guided)
  for (int i=0; i< natom; i++) NumNei[i] = 0;

  for (int k=0; k < natom; k++){
    int ik = map[atom->type[k]];
    double Epk = 0.;
    #pragma omp parallel for default(shared) schedule(guided) reduction(+:Epk)
    for (int kp=k; kp < natom; kp++){
      int ikp = map[atom->type[kp]];
      for (int kx = -Lx; kx <= Lx; kx++)
      for (int ky = -Ly; ky <= Ly; ky++)
      for (int kz = -Lz; kz <= Lz; kz++){
        if ( (k == kp) && (kx==0 && ky==0 && kz==0) ) continue;
        double Rklkp[3];
        Rklkp[0] = atom->s[kp][0] - atom->s[k][0] + double(kx);
        Rklkp[1] = atom->s[kp][1] - atom->s[k][1] + double(ky);
        Rklkp[2] = atom->s[kp][2] - atom->s[k][2] + double(kz);

        double r2 = atom->veclen2(Rklkp);
        if (r2 >= eam->rcutsq) continue;

        double r = sqrt(r2);
        #pragma omp atomic
        den[k] += eam->Rho(r, ikp, ik);
        if (k != kp){
          den[kp] += eam->Rho(r, ik, ikp);
          Epk += eam->Phi(r, ik, ikp);
        } else {
          Epk += 0.5*eam->Phi(r, ikp, ik);
        }

        #pragma omp critical
        { // only one thread is allowed to modify this region at the same time
          if (++NumNei[k] > neimax){
            neimax += 12;
            memory->grow(NeiList,neimax,natom,"DYNMAT_setup:NeiList");
            memory->grow(Bonds,neimax,natom,4,"DYNMAT_setup:Bonds");
          }
  
          int nnow = NumNei[k] -1;
          NeiList[nnow][k] = kp;
          for (int idim=0; idim<3; idim++) Bonds[nnow][k][idim] = Rklkp[idim];
          Bonds[nnow][k][3] = r;
        }  
        if (k != kp){
          if ( ++NumNei[kp]> neimax){
            neimax += 12;
            memory->grow(NeiList,neimax,natom,"DYNMAT_setup:NeiList");
            memory->grow(Bonds,neimax,natom,4,"DYNMAT_setup:Bonds");
          }
          int nnow = NumNei[kp] -1;
          NeiList[nnow][kp] = k;
          for (int idim=0; idim<3; idim++) Bonds[nnow][kp][idim] = -Rklkp[idim];
          Bonds[nnow][kp][3] = r;
        }
      }
    }
    Ep += Epk;
  }
  // get the embedded energy and its derivatives
  double Eci = 0.;
#pragma omp parallel for default(shared) schedule(guided) reduction(+:Eci)
  for (int i=0; i<natom; i++){
    int ip  = map[atom->type[i]];
    Eci += eam->F(den[i], ip);

    Fp[i]  = eam->Fp(den[i], ip);
    Fpp[i] = eam->Fpp(den[i], ip);
  }
  Ec = Eci;
  Et = Ec + Ep;
  // display lattice energy
  printf("\n");for (int i=0; i<60;i++) printf("="); printf("\n");
  printf("System name : %s\n", atom->title);
  printf("System size : %d atoms with %d types\n",natom, atom->ntype);
  printf("Elements are:");for (int i=0; i<atom->ntype;i++) printf(" %s -> %s,",atom->elements[i], eam->elements[map[i]]);
  printf("\nCohesive eng: %lg eV\n", Ec);
  printf("Pairwise eng: %lg eV\n", Ep);
  printf("Total energy: %lg eV ==> %lg eV/atom\n", Et, Et/double(natom));
  for (int i=0; i<60;i++) printf("="); printf("\n");
}

/*------------------------------------------------------------------------------
 * Method to calculate the dynamical matrix at q
 *----------------------------------------------------------------------------*/
void DYNMAT::computeDM(double *q)
{
  const double tpi = 8.*atan(1.);
  for (int i=0; i<ndim; i++){
    for (int j=0; j<ndim; j++){dm[i][j].REALPART = 0.; dm[i][j].IMAGPART=0.;}
  }

  // now to compute the dynamical matrix
  #pragma omp parallel for default(shared) schedule(guided)
  for (int k=0; k<natom; k++){
    double Dkk[3][3];
    for (int i=0; i<3; i++){Dkk[i][0] = Dkk[i][1] = Dkk[i][2] = 0.;}
    int ik = map[atom->type[k]];

    for (int m=0; m<NumNei[k]; m++){
      int kp = NeiList[m][k];
      int ikp = map[atom->type[kp]];
      
      double Rkkp[3], Xkkp[3];
      double r, r_inv, r2_inv, r3_inv;
      for (int idim=0; idim<3; idim++) Rkkp[idim] = Bonds[m][k][idim];
      r = Bonds[m][k][3];
      r_inv = 1./r; r2_inv = r_inv*r_inv; r3_inv = r_inv * r2_inv;

      atom->dir2car(Rkkp, Xkkp);

      double fpk, fpkp, fppk, fppkp, phip, phipp;
      fpk   = eam->Rhop(r, ikp, ik);
      fpkp  = eam->Rhop(r, ik, ikp);
      fppk  = eam->Rhopp(r, ikp, ik);
      fppkp = eam->Rhopp(r, ik, ikp);

      phip  = eam->Phip(r, ikp, ik);
      phipp = eam->Phipp(r, ikp, ik);

      double qr = tpi*(q[0]*Rkkp[0]+q[1]*Rkkp[1]+q[2]*Rkkp[2]);
      double rmass = 1./sqrt(eam->mass[ik]*eam->mass[ikp]);
      COMPLEX iqr;

      iqr.REALPART = cos(qr)*rmass;
      iqr.IMAGPART = sin(qr)*rmass;

      double Term_ab = -(Fpp[k]*fpk *fpk + Fp[k]*fppk + Fpp[kp]*fpkp*fpkp + Fp[kp]*fppkp) * r2_inv
                       +(Fp[k]*fpk + Fp[kp]*fpkp) * r3_inv - phipp * r2_inv + phip * r3_inv;
      double Term_aa = -(Fp[k]*fpk + Fp[kp]*fpkp + phip)*r_inv;

      for (int a=0; a<3; a++)
      for (int b=0; b<3; b++){
        double RaRb = Xkkp[a]*Xkkp[b];
        double Dkakpb = Term_ab * RaRb;
        if (a == b) Dkakpb += Term_aa;

        int i = a + 3*k, j = b + 3*kp;
        dm[i][j].REALPART += Dkakpb * iqr.REALPART;
        dm[i][j].IMAGPART += Dkakpb * iqr.IMAGPART;

        Dkk[a][b] -= Dkakpb * rmass;
      }
    }
    int ii = 3*k;
    for (int i=0; i<3; i++){
      for (int j=0; j<3; j++) dm[ii+i][ii+j].REALPART += Dkk[i][j];
    }
  }
return;
}

/*------------------------------------------------------------------------------
 * Method to calculate the eigen values and eigen vectors of the dynamical matrix
 *----------------------------------------------------------------------------*/
int DYNMAT::computeEigen(int flag)
{
  char jobz, uplo;
  const double eva2thz = sqrt(1.60217733*6.022142e3); // works if the original units are eV and Angstrom

  INTEGER n, lda, lwork, lrwork, *iwork, liwork, info;
  COMPLEX *work;
  DOUBLEREAL *w = &egval[0], *rwork;

  if (flag) jobz = 'V';
  else jobz = 'N';
  uplo = 'U';
  n     = ndim;
  lwork = (n + 2)*n;
  lda    = n;
  lrwork = 1 + (5+n+n)*n;
  liwork = 3 + 5*n;
  memory->create(work,  lwork,"computeEigen:work");
  memory->create(rwork, lrwork,"computeEigen:rwork");
  memory->create(iwork, liwork,"computeEigen:iwork");

#ifdef OMP
#ifdef MKL
  int npmax = omp_get_max_threads();
  omp_set_num_threads(npmax);
#endif
#endif

  zheevd_(&jobz, &uplo, &n, dm[0], &lda, w, work, &lwork, rwork, &lrwork, iwork, &liwork, &info);

  // to get w instead of w^2; and convert the unit of w
  for (int i=0; i<n; i++){
    if (w[i]>= 0.) w[i] = sqrt(w[i]);
    else w[i] = -sqrt(-w[i]);

    w[i] *= eva2thz;
  }

  memory->destroy(work);
  memory->destroy(rwork);
  memory->destroy(iwork);

return info;
}

/*------------------------------------------------------------------------------
 * Private method to choose the EAM parameterizations
 *----------------------------------------------------------------------------*/
void DYNMAT::selectEAM(void)
{
  // Main menu to select EAM parameterization
  int pottype = 1;
  char str[MAXLINE], *ptr;

  printf("\n"); for (int i=0; i<60; i++) printf("="); printf("\n");
  printf("Please select the parameterization of the EAM potential:\n");
  printf(" 1. Setfl format;\n");
  printf(" 2. Funcfl format;\n");
  printf(" 3. EAM/Finnis-Sinclair;\n");
  printf("Your choice [1]: ");
  fgets(str,MAXLINE,stdin);
  ptr = strtok(str, " \t\n\r\f");
  if (ptr) pottype = atoi(ptr);
  
  if (pottype == 1) eam = new SETFL;
  else if (pottype == 2) eam = new FUNCFL;
  else if (pottype == 3) eam = new EAMFS;
  else{
    printf("\nError: Wrong type of EAM selected!\n");
    exit(1);
  }

  if (eam->initialized == 0){
    printf("\nError encountered while initializing EAM!\n");
    exit(2);
  }
  eam->check_eam();
}

/*------------------------------------------------------------------------------
 * Private method to calculate the LDOS by using the Green function method
 *
 * nt  : on the order of 100 is usually enough, for large system, should be 
 *       increased slightly;
 * eps : on the order of 10~100, small for large system.
 *----------------------------------------------------------------------------*/
void DYNMAT::GreenLDOS()
{
  double q0[3], wmin, wmax, eps;
  int sysdim = 3, nt = ndim/10, ndos = 201;
  char str[MAXLINE];
  const double eva2thz =1.60217733*6.022142e3; // assuming original units are eV and Angstrom
 
  printf("\nNow to evaluate the Hessian at Gamma-point ...");
  double **hessian = (double **) dm;
  // to get the dynamical matrix at gamma point
  q0[0] = q0[1] = q0[2] = 0.;
  computeDM(&q0[0]);

  for (int idim=0; idim<ndim; idim++){
    for (int jdim=0; jdim<ndim; jdim++) hessian[idim][jdim] = dm[idim][jdim].REALPART * eva2thz;
  }
  printf("Done!\n");

  wmin = 0.;
  wmax = 10.;
  eps  = 12.;

  int nlocal = 0, *locals;
  printf("\nThe total number of atoms in your system is %d.\n", natom);
  printf("Please input the atom index to get local PDOS [1]: ");
  int nwd = atom->count_words(fgets(str,MAXLINE,stdin));
  if (nwd < 1){nwd = 1; strcpy(str,"1");}

  char *ptr = strtok(str, " \t\n\r\f");
  char id4ldos[MAXLINE];

  if (strcmp(ptr, "=") == 0){
    strcpy(id4ldos, "=");
    if (nwd > 1){
      memory->create(locals, nwd-1, "GreenLDOS:locals");

      ptr = strtok(NULL," \t\n\r\f");
      while (ptr){
        int id = atoi(ptr);
        if (id > 0 && id <= natom){
          int flag = 1;
          for (int i=0; i<nlocal; i++){
            if (id == locals[i]){ flag = 0; break; }
          }
          if (flag){
            locals[nlocal++] = id;
            strcat(id4ldos, " "); strcat(id4ldos, ptr); 
          }
        }
  
        ptr = strtok(NULL," \t\n\r\f");
      }
    }
  } else if (strcmp(ptr, ">") == 0){
    if (nwd > 1){
      int nlo = atoi(strtok(NULL," \t\n\r\f"));
      nlocal = natom-nlo;

      if (nlocal > 0 && nlocal <= natom){
        memory->create(locals, nlocal, "GreenLDOS:locals");
        for (int i=0; i<nlocal; i++) locals[i] = i+nlo+1;

        sprintf(id4ldos,"> %d", nlo);
      } else nlocal = 0;
    }
  } else if (strcmp(ptr, ">=") == 0){
    if (nwd > 1){
      int nlo = atoi(strtok(NULL," \t\n\r\f"));
      nlocal = natom-nlo+1;
      if (nlocal > 0 && nlocal <= natom){
        memory->create(locals, nlocal, "GreenLDOS:locals");
        for (int i=0; i<nlocal; i++) locals[i] = i+nlo;

        sprintf(id4ldos,">= %d", nlo);
      } else nlocal = 0;
    }
  } else if (strcmp(ptr, "<") == 0){
    if (nwd > 1){
      int nhi = atoi(strtok(NULL," \t\n\r\f"));
      nlocal = nhi-1;
      if (nlocal > 0 && nlocal <= natom){
        memory->create(locals, nlocal, "GreenLDOS:locals");
        for (int i=0; i<nlocal; i++) locals[i] = i;

        sprintf(id4ldos,"< %d", nhi);
      } else nlocal = 0;
    }
  } else if (strcmp(ptr, "<=") == 0){
    if (nwd > 1){
      int nhi = atoi(strtok(NULL," \t\n\r\f"));
      nlocal = nhi;
      if (nlocal > 0 && nlocal <= natom){
        memory->create(locals, nlocal, "GreenLDOS:locals");
        for (int i=0; i<nlocal; i++) locals[i] = i;

        sprintf(id4ldos,"<= %d", nhi);
      } else nlocal = 0;
    }
  } else if (strcmp(ptr, "<>") == 0){
    if (nwd > 2){
      int nlo = atoi(strtok(NULL," \t\n\r\f"));
      int nhi = atoi(strtok(NULL," \t\n\r\f"));

      if (nlo > 0 && nhi <= natom && nhi >= nlo){
        nlocal = nhi - nlo + 1;
        memory->create(locals, nlocal, "GreenLDOS:locals");

        for (int i=0; i<nlocal; i++) locals[i] = i+nlo;

        sprintf(id4ldos,"[%d, %d]", nlo, nhi);
      }
    }
  } else if (strcmp(ptr, "><") == 0){
    if (nwd > 2){
      int nlo = atoi(strtok(NULL," \t\n\r\f"));
      int nhi = atoi(strtok(NULL," \t\n\r\f"));

      if (nlo > 0 && nhi <= natom && nhi >= nlo){
        nlocal = nlo + natom - nhi + 1;
        memory->create(locals, nlocal, "GreenLDOS:locals");

        for (int i=0;     i<= nlo;  i++) locals[i] = i;
        for (int i=nlo+1; i<nlocal; i++) locals[i] = i-nlo-1+nhi;

        sprintf(id4ldos,"[1, %d] & [%d, %d]", nlo, nhi, natom);
      }
    }
  } else {
    nlocal = nwd;
    if (nlocal > 0){
      memory->create(locals, nlocal, "GreenLDOS:locals");
      nlocal = 0;

      strcpy(id4ldos, "=");
      while (ptr){
        int id = atoi(ptr);
        if (id > 0 && id <= natom){
          int flag = 1;
          for (int i=0; i<nlocal; i++){
            if (id == locals[i]){ flag = 0; break; }
          }
          if (flag){
            locals[nlocal++] = id;
            strcat(id4ldos, " "); strcat(id4ldos, ptr);
          }
        }

        ptr = strtok(NULL," \t\n\r\f");
      }
    }
  }
  if (nlocal < 1){
    hessian = NULL;
    return;
  }
  printf("Local phonon DOS for %d atoms with id: %s will be computed!", nlocal, id4ldos);

  printf("\nPlease input the maximum iteration during Lanczos [%d]: ", nt);
  if (atom->count_words(fgets(str,MAXLINE,stdin)) > 0) nt = atoi(strtok(str," \t\n\r\f"));
  if (nt < 4){
    printf("\nError: too few Lanczos steps!\n");
    hessian = NULL;
    return;
  }

  printf("Please input the number of points for LDOS [%d]: ", ndos);
  if (atom->count_words(fgets(str,MAXLINE,stdin)) > 0) ndos = atoi(strtok(str," \t\n\r\f"));
  if (ndos < 10){
    printf("\nError: too few DOS points!\n");
    hessian = NULL;
    return;
  }
  ndos += (ndos+1)%2;

  printf("Please input the frequency range to evaluate LDOS [%lg %lg]: ", wmin, wmax);
  if (atom->count_words(fgets(str,MAXLINE,stdin)) > 1){
    wmin = atof(strtok(str," \t\n\r\f"));
    wmax = atof(strtok(NULL," \t\n\r\f"));
  }
  if (wmax < wmin || wmax < 0.){
    printf("\nError: incorrect frequency range!\n");
    hessian = NULL;
    return;
  }

  printf("Please input the value of epsilon in real space Green function [%lg]: ", eps);
  if (atom->count_words(fgets(str,MAXLINE,stdin)) > 0) eps = atof(strtok(str," \t\n\r\f"));
  if (eps <= 0.){
    printf("\nError: wrong epsilon, expected a positive number!\n");
    hessian = NULL;
    return;
  }

  Timer *timer = new Timer();
#ifdef OMP
  int npmax = omp_get_max_threads();
#endif
  #pragma omp parallel for default(shared) schedule(guided) num_threads(MIN(nlocal,npmax))
  for (int ilocal = 0; ilocal<nlocal; ilocal++){
    int iatom = locals[ilocal];
#ifndef OMP
    printf("Now to compute PLDOS for atom %d... ", iatom);
#else
    Timer *timer_in = new Timer();
#endif

    // now to compute the LDOS by using class Green
    Green *green = new Green(natom, sysdim, nt, wmin, wmax, ndos, eps, hessian, iatom);
    delete green;

#ifndef OMP
    printf("Done! Wall time used: %g/%g seconds.\n", timer->sincelast(), timer->up2now());
#else
    int me = omp_get_thread_num();
    timer_in->stop();
    printf("PLDOS of atom %d computed by %d! Wall time used: %g seconds.\n", iatom, me, timer_in->wall_time() );
    delete timer_in;
#endif
    fflush(stdout);
  }
  timer->stop();
  printf("Total CPU time used: %g seconds; wall time: %g seconds.\n", timer->cpu_time(), timer->wall_time());
  delete timer;

  hessian = NULL;

return;
}
