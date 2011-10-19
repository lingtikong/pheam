#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "timer.h"

#include "dynmat.h"

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
  map = memory->create(map,atom->ntype,"DYNMAT:map");
  checkmap();

  // allocate variables
  natom = atom->natom;
  ndim = natom*3;
  neimax = 40;

  den = memory->create(den, natom, "DYNMAT:den");
  Fp  = memory->create(Fp , natom, "DYNMAT:Fp");
  Fpp = memory->create(Fpp, natom, "DYNMAT:Fpp");
  NumNei = memory->create(NumNei, natom, "DYNMAT:NumNei");

  NeiList = memory->create(NeiList,neimax, natom,"DYNMAT_DYNMAT:NeiList");
  Bonds   = memory->create(Bonds,neimax,natom,4,"DYNMAT_DYNMAT:Bonds");

  // ask for the extension of lattice in three dimensions
  // assuming orthogonal lattice
  Lx = int(eam->rcut/atom->axis[0][0]+1.0); if (atom->pbc[0] == 0) Lx=0;
  Ly = int(eam->rcut/atom->axis[1][1]+1.0); if (atom->pbc[1] == 0) Ly=0;
  Lz = int(eam->rcut/atom->axis[2][2]+1.0); if (atom->pbc[2] == 0) Lz=0;
  char str[MAXLINE];
  printf("\nPlease input the # of extension of cell in three dimensions[%d %d %d]: ", Lx,Ly,Lz);
  if (eam->count_words(gets(str)) >= 3) sscanf(str,"%d %d %d", &Lx, &Ly, &Lz);
  if (atom->pbc[0] == 0) Lx=0;
  if (atom->pbc[1] == 0) Ly=0;
  if (atom->pbc[2] == 0) Lz=0;
  printf("The # of extension in three dimensions accepted: %d x %d x %d.", Lx,Ly,Lz);

  // now to get the energy, electronic density and embedded terms of the system
  setup();

  // allocate memory for dynamical matrix and its eigen values/vectors
  egval   = memory->create(egval,ndim,"DYNMAT_DYNMAT:egval");
  dm      = memory->create(dm,ndim, ndim,"DYNMAT_DYNMAT:dm");

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
  for (int i=0; i<atom->ntype; i++){
    int ip = eam->index(atom->elements[i]);
    if (ip < 0){
      printf("\nError: Cannot find element %s in EAM; available elements in EAM are:", atom->elements[i]);
      for (int j=0; j<eam->ntype; j++) printf(" %s", eam->elements[j]); printf("\n");
      while (ip<0){
        printf("Please input the EAM name for %s [%s]: ", atom->elements[i], eam->elements[0]);
        char ename[10];
        if (strlen(gets(ename)) < 1) strcpy(ename, eam->elements[0]);
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
  for (int i=0; i< natom; i++) den[i] = 0.;
  for (int i=0; i< natom; i++) NumNei[i] = 0;

  for (int k=0; k < natom; k++){
    int ik = map[atom->type[k]];
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

        double r = sqrt(r2); //, rhokpk = eam->Rho(r, ikp, ik);;
        //den[k] += rhokpk;
        //if (k == 512 || kp == 512) printf("k=%d kp=%d, r=%g rho=%g\n", k, kp, r, rhokpk);
        den[k] += eam->Rho(r, ikp, ik);
        if (k != kp){den[kp] += eam->Rho(r, ik, ikp); Ep += eam->Phi(r, ik, ikp);}
        else Ep += 0.5*eam->Phi(r, ikp, ik);

        if (++NumNei[k] > neimax){
          neimax += 12;
          NeiList = memory->grow(NeiList,neimax,natom,"DYNMAT_setup:NeiList");
          Bonds   = memory->grow(Bonds,neimax,natom,4,"DYNMAT_setup:Bonds");
        }
        int nnow = NumNei[k] -1;
        NeiList[nnow][k] = kp;
        for (int idim=0; idim<3; idim++) Bonds[nnow][k][idim] = Rklkp[idim];
        Bonds[nnow][k][3] = r;
          
        if (k != kp){
          if ( ++NumNei[kp]> neimax){
            neimax += 12;
            NeiList = memory->grow(NeiList,neimax,natom,"DYNMAT_setup:NeiList");
            Bonds   = memory->grow(Bonds,neimax,natom,4,"DYNMAT_setup:Bonds");
          }
          int nnow = NumNei[kp] -1;
          NeiList[nnow][kp] = k;
          for (int idim=0; idim<3; idim++) Bonds[nnow][kp][idim] = -Rklkp[idim];
          Bonds[nnow][kp][3] = r;
        }
      }
    }
  }
  // get the embedded energy and its derivatives
  for (int i=0; i<natom; i++){
    int ip  = map[atom->type[i]];
    //printf("i=%d, ip=%d, rho=%g\n", i, ip, den[i]);
    Ec += eam->F(den[i], ip);

    Fp[i]  = eam->Fp(den[i], ip);
    Fpp[i] = eam->Fpp(den[i], ip);
  }
  Et = Ec + Ep;
  // display lattice energy
  printf("\n");for (int i=0; i<60;i++) printf("="); printf("\n");
  printf("System name : %s\n", atom->title);
  printf("System size : %d atoms with %d types\n",natom, atom->ntype);
  printf("Elements are:");for (int i=0; i<atom->ntype;i++) printf(" %s",atom->elements[i]);
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
    for (int j=0; j<ndim; j++){dm[i][j].r = 0.; dm[i][j].i=0.;}
  }

  // now to compute the dynamical matrix
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
      doublecomplex iqr;

      iqr.r = cos(qr)*rmass;
      iqr.i = sin(qr)*rmass;

      double Term_ab = -(Fpp[k]*fpk *fpk + Fp[k]*fppk + Fpp[kp]*fpkp*fpkp + Fp[kp]*fppkp) * r2_inv
                       +(Fp[k]*fpk + Fp[kp]*fpkp) * r3_inv - phipp * r2_inv + phip * r3_inv;
      double Term_aa = -(Fp[k]*fpk + Fp[kp]*fpkp + phip)*r_inv;

      for (int a=0; a<3; a++)
      for (int b=0; b<3; b++){
        double RaRb = Xkkp[a]*Xkkp[b];
        double Dkakpb = Term_ab * RaRb;
        if (a == b) Dkakpb += Term_aa;

        int i = a + 3*k, j = b + 3*kp;
        dm[i][j].r += Dkakpb * iqr.r;
        dm[i][j].i += Dkakpb * iqr.i;

        Dkk[a][b] -= Dkakpb * rmass;
      }
    }
    int ii = 3*k;
    for (int i=0; i<3; i++){
      for (int j=0; j<3; j++) dm[ii+i][ii+j].r += Dkk[i][j];
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
  integer n, lda, lwork, lrwork, *iwork, liwork, info;
  doublecomplex *work;
  doublereal *w = &egval[0], *rwork;

  if (flag) jobz = 'V';
  else jobz = 'N';
  uplo = 'U';
  n     = ndim;
  lwork = (n + 2)*n;
  lda    = n;
  lrwork = 1 + (5+n+n)*n;
  liwork = 3 + 5*n;
  work  = memory->create(work, lwork,"computeEigen:work");
  rwork = memory->create(rwork, lrwork,"computeEigen:rwork");
  iwork = memory->create(iwork, liwork,"computeEigen:iwork");

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
  char str[MAXLINE];

  printf("\n"); for (int i=0; i<60; i++) printf("="); printf("\n");
  printf("Please select the parameterization of the EAM potential:\n");
  printf(" 1. Setfl format;\n");
  printf(" 2. Funcfl format;\n");
  printf(" 3. EAM/Finnis-Sinclair;\n");
  printf("Your choice [1]: ");
  if (strlen(gets(str)) > 0) sscanf(str,"%d", &pottype);
  
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
  double **hessian, q0[3], wmin, wmax, eps;
  int sysdim = 3, nt = ndim*0.1, iatom = 1, ndos = 201;
  char str[MAXLINE];
  const double eva2thz =1.60217733*6.022142e3; // assuming original units are eV and Angstrom
 
  printf("\nNow to evaluate the Hessian at Gamma-point ...");
  hessian = memory->create(hessian,ndim,ndim,"DYNMAT_GreenLDOS:hessian");
  // to get the frequency range at gamma point, as well as the dynamical matrix
  q0[0] = q0[1] = q0[2] = 0.;
  computeDM(&q0[0]);

  for (int idim=0; idim<ndim; idim++){
    for (int jdim=0; jdim<ndim; jdim++) hessian[idim][jdim] = dm[idim][jdim].r * eva2thz;
  }
  printf("Done!\n");

  wmin = 0.;
  wmax = 10.;
  eps  = 12.;

  Timer *time = new Timer();

  while ( 1 ) {
    printf("\nThe total number of atoms in your system is %d.\n", natom);
    printf("Please input the atom index to get local PDOS [%d]: ", iatom);
    if (atom->count_words(gets(str)) > 0) iatom = atoi(strtok(str," \t\n\r\f"));
    if (iatom < 1 || iatom > natom) break;

    printf("Please input the maximum iteration during Lanczos [%d]: ", nt);
    if (strlen(gets(str)) > 0) nt = atoi(strtok(str," \t\n\r\f"));
    if (nt < 4) continue;

    printf("Please input the number of points for LDOS [%d]: ", ndos);
    if (strlen(gets(str)) > 0) ndos = atoi(strtok(str," \t\n\r\f"));
    if (ndos < 10) continue;
    ndos += (ndos+1)%2;

    printf("Please input the frequency range to evaluate LDOS [%lg %lg]: ", wmin, wmax);
    if (atom->count_words(gets(str)) > 1){wmin = atof(strtok(str," \t\n\r\f")); wmax = atof(strtok(NULL," \t\n\r\f"));}
    if (wmax < wmin || wmax < 0.) continue;

    printf("Please input the value of epsilon in real space Green function [%lg]: ", eps);
    if (strlen(gets(str)) > 0) eps = atof(strtok(str," \t\n\r\f"));
    if (eps <= 0.) continue;

    time->start();
    // now to compute the LDOS by using class Green
    green = new Green(natom, sysdim, nt, wmin, wmax, ndos, eps, hessian, iatom);
    delete green;
    time->stop(); time->print();
  }

  memory->destroy(hessian);
return;
}
