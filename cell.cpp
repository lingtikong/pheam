#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "cell.h"

#define MAXLINE 512

//using namespace std;

/*------------------------------------------------------------------------------
 * Constructor is used to read in the atomic position file, which is assumed to
 * be in the xyz format.
 *----------------------------------------------------------------------------*/
CELL::CELL()
{
  // create memory and elements array
  memory = new Memory();
  elements = memory->create(elements,10,10,"CELL_CELL:elements");

  // ask for file name
  char str[MAXLINE], *ptr;
  while (1){
    printf("\nPlease input the xyz file name: ");
    fgets(str, MAXLINE, stdin);
    ptr = strtok(str, " \t\n\r\f");
    if (ptr) break;
  }

  int n = strlen(ptr) + 1;
  char *fname = new char[n];
  strcpy(fname, ptr);

  // open atomic position file
  FILE *fp = fopen(fname, "r");
  if (fp == NULL){printf("\nError: file %s not found!\n", fname); exit(1);}

  // read atomic position file 
  fgets(str,MAXLINE,fp);
  sscanf(str, "%d", &natom);
  fgets(str,MAXLINE,fp);
  n = strlen(str) + 1;
  title = new char[n];
  strcpy(title, str);
  if (natom < 1){
    printf("\nError: wrong number of atoms (%d) read from file: %s\n", natom, fname);
    exit(2);
  }
  x = memory->create(x,natom,3,"CELL_CELL:x");
  s = memory->create(s,natom,3,"CELL_CELL:s");
  axis    = memory->create(axis,3,3,"CELL_CELL:axis");
  invaxis = memory->create(invaxis,3,3,"CELL_CELL:invaxis");
  type = memory->create(type, natom,"CELL_CELL:type");

  char *name, *sdum;
  int flag[3];
  ntype = 0;
  flag[0] = flag[1] = flag[2] = 0;
  for (int i=0; i< natom; i++){
    int dir;
    if (fgets(str,MAXLINE,fp) == NULL){
      printf("\nError while reading atomic positions from file: %s\n", fname);
      exit(3);
    }
    int nr = count_words(str);
    if (nr >= 4){
      name = strtok(str," \t\n\r\f");
      type[i] = index(name);
      x[i][0] = atof(strtok(NULL," \t\n\r\f"));
      x[i][1] = atof(strtok(NULL," \t\n\r\f"));
      x[i][2] = atof(strtok(NULL," \t\n\r\f"));
      if (nr >= 9){
        sdum = strtok(NULL," \t\n\r\f");
        dir  = atoi(strtok(NULL," \t\n\r\f"))-1;

        if (dir >=0 && dir <3){
          axis[dir][0] = atof(strtok(NULL," \t\n\r\f"));
          axis[dir][1] = atof(strtok(NULL," \t\n\r\f"));
          axis[dir][2] = atof(strtok(NULL," \t\n\r\f"));
          flag[dir] = 1;
        }
      }
    } else {
      printf("\nError while reading the %d-th atom from file: %s\n", i+1, fname);
      exit(3);
    }
  }
  fclose(fp);

  // ask for lattice info if not read correctly from file
  if (flag[0]==0 || flag[1]==0 || flag[2]==0){
    printf("\nLattice info read from %s is insufficient, please input them properly.\n", fname);
    for (int i=0; i<3; i++){
      do printf("Please input the vector A%d: ", i+1);
      while (count_words(fgets(str,MAXLINE,stdin)) < 3);
      axis[i][0] = atof(strtok(str," \t\n\r\f"));
      axis[i][1] = atof(strtok(NULL," \t\n\r\f"));
      axis[i][2] = atof(strtok(NULL," \t\n\r\f"));
    }
  }
  // ask for boundary condition
  printf("Please indicate if PBC is applied in each direction (1:yes; 0:no)[1 1 1]: ");
  if (count_words(fgets(str,MAXLINE,stdin)) >= 3){
    pbc[0] = atoi(strtok(str, " \t\n\r\f"));
    pbc[1] = atoi(strtok(NULL," \t\n\r\f"));
    pbc[2] = atoi(strtok(NULL," \t\n\r\f"));
  } else pbc[0] = pbc[1] = pbc[2] = 1;

  for (int i=0; i<3; i++) if (pbc[i] != 0) pbc[i] = 1;

  // get fractional coordinate as well.
  GaussJordan(3,axis[0],invaxis[0]);
  car2dir();
  delete []fname;

return;
}

/*------------------------------------------------------------------------------
 * Deconstructor is used to free all allocated memories.
 *----------------------------------------------------------------------------*/
CELL::~CELL()
{
  memory->destroy(elements);
  memory->destroy(x);
  memory->destroy(s);
  memory->destroy(axis);
  memory->destroy(invaxis);
  memory->destroy(type);

  delete []title;
  delete memory;
}

/*------------------------------------------------------------------------------
 * Method to convert cartesian coordinate into fractional
 *----------------------------------------------------------------------------*/
void CELL::car2dir()
{
  for (int i=0; i<natom; i++){
    for (int idim=0; idim<3; idim++){
      s[i][idim] = 0.;
      for (int jdim=0; jdim<3; jdim++) s[i][idim] += x[i][jdim]*invaxis[jdim][idim];

      while (s[i][idim] >= 1.) s[i][idim] -= 1.;
      while (s[i][idim] <  0.) s[i][idim] += 1.;
    }
  }

return;
}

/*------------------------------------------------------------------------------
 * Method to convert fractional coordinate into cartesian
 *----------------------------------------------------------------------------*/
void CELL::dir2car()
{
  for (int i=0; i<natom; i++){
    for (int idim=0; idim<3; idim++){
      x[i][idim] = 0.;
      for (int jdim=0; jdim<3; jdim++) x[i][idim] += s[i][jdim]*axis[jdim][idim];
    }
  }

return;
}

/*------------------------------------------------------------------------------
 * Method to convert a fractional vector into a cartesian one
 *----------------------------------------------------------------------------*/
void CELL::dir2car(double *frac, double *cart)
{
  for (int idim=0; idim<3; idim++){
    cart[idim] = 0.;
    for (int jdim=0; jdim<3; jdim++) cart[idim] += frac[jdim]*axis[jdim][idim];
  }
return;
}

/*------------------------------------------------------------------------------
 * Method to calculate the length of a fractional vector
 *----------------------------------------------------------------------------*/
double CELL::veclen(double *ss)
{
  double xx[3];
  for (int idim=0; idim<3; idim++){
    xx[idim] = 0.;
    for (int jdim=0; jdim<3; jdim++) xx[idim] += ss[jdim]*axis[jdim][idim];
  }
  return sqrt(xx[0]*xx[0]+xx[1]*xx[1]+xx[2]*xx[2]);
}

/*------------------------------------------------------------------------------
 * Method to calculate the square of the length of a fractional vector
 *----------------------------------------------------------------------------*/
double CELL::veclen2(double *ss)
{
  double xx[3];
  for (int idim=0; idim<3; idim++){
    xx[idim] = 0.;
    for (int jdim=0; jdim<3; jdim++) xx[idim] += ss[jdim]*axis[jdim][idim];
  }
  return xx[0]*xx[0]+xx[1]*xx[1]+xx[2]*xx[2];
}

/*------------------------------------------------------------------------------
 * Method to calculate the length of a cartesian vector
 *----------------------------------------------------------------------------*/
double CELL::VecLen(double *xx)
{
  return sqrt(xx[0]*xx[0]+xx[1]*xx[1]+xx[2]*xx[2]);
}

/*------------------------------------------------------------------------------
 * Method to calculate the square of the length of a cartesian vector
 *----------------------------------------------------------------------------*/
double CELL::VecLen2(double *xx)
{
  return xx[0]*xx[0]+xx[1]*xx[1]+xx[2]*xx[2];
}

/*------------------------------------------------------------------------------
 * Method to calculate the cosine of the angle beteen two vectors
 *----------------------------------------------------------------------------*/
double CELL::VecAng(double *xx, double *yy)
{
  double lx = VecLen(xx), ly = VecLen(yy);
  return (xx[0]*yy[0]+xx[1]*yy[1]+xx[2]*yy[2])/(lx*ly);
}

/*------------------------------------------------------------------------------
 * Private method to get the index of the atomic type
 *----------------------------------------------------------------------------*/
int CELL::index(const char *name)
{
  int ip=-1;
  for (int i=0; i<ntype; i++){
    if (strcmp(name, elements[i]) == 0){ip = i; break;}
  }
  if (ip < 0){
    ip = ntype;
    strcpy(elements[ntype++], name);
  }
  return ip;
}

/*------------------------------------------------------------------------------
 * Method to display the lattice info
 *----------------------------------------------------------------------------*/
void CELL::display( )
{
  printf("\n");for (int i=0; i<60; i++) printf("=");printf("\n");
  printf("Total number of atomic types in cell : %d\n", ntype);
  printf("The names of these elements  are     :");
  for (int i=0; i<ntype; i++) printf(" %s", elements[i]);
  printf("\nTotal number of atoms in unit cell   : %d\n", natom);
  printf("The basis vectors of the cell:\n");
  for (int i=0; i<3; i++) printf("A%d: %lg %lg %lg\n", i+1, axis[i][0], axis[i][1], axis[i][2]);
  for (int i=0; i<60; i++) printf("=");printf("\n");

return;
}

/*------------------------------------------------------------------------------
 * Method to count # of words in a string, without destroying the string
 *----------------------------------------------------------------------------*/
int CELL::count_words(const char *line)
{
  int n = strlen(line) + 1;
  char *copy = memory->create(copy, n,"copy");
  strcpy(copy,line);

  char *ptr;
  if (ptr = strchr(copy,'#')) *ptr = '\0';

  if (strtok(copy," \t\n\r\f") == NULL) {
    memory->destroy(copy);
    return 0;
  }
  n = 1;
  while (strtok(NULL," \t\n\r\f")) n++;

  memory->destroy(copy);
  return n;
}

/*------------------------------------------------------------------------------
 * Private method to do matrix inversion
 *----------------------------------------------------------------------------*/
void CELL::GaussJordan(const int n, const double *MatA, double *Mat)
{
  int i,icol,irow,j,k,l,ll,idr,idc;
  int indxc[n],indxr[n],ipiv[n];
  double big, dum, pivinv;

  for (int i=0; i<n*n; i++) Mat[i] = MatA[i];

  for (i=0; i<n; i++) ipiv[i] = 0;
  for (i=0; i<n; i++){
    big = 0.;
    for (j=0; j<n; j++){
      if (ipiv[j] != 1){
        for (k=0; k<n; k++){
          if (ipiv[k] == 0){
            idr = j*n+k;
            if (fabs(Mat[idr]) >= big){
              big  = fabs(Mat[idr]);
              irow = j;
              icol = k;
            }
          }else if (ipiv[k] >1){
            printf("\nError: Singular matrix in double GaussJordan!\n");
          }
        }
      }
    }
    ipiv[icol] += 1;
    if (irow != icol){
      for (l=0; l<n; l++){
        idr  = irow*n+l;
        idc  = icol*n+l;
        dum  = Mat[idr];
        Mat[idr] = Mat[idc];
        Mat[idc] = dum;
      }
    }
    indxr[i] = irow;
    indxc[i] = icol;
    idr = icol*n+icol;
    if (Mat[idr] == 0.) printf("\nError: Singular matrix in double GaussJordan!\n");
    
    pivinv = 1./ Mat[idr];
    Mat[idr] = 1.;
    idr = icol*n;
    for (l=0; l<n; l++) Mat[idr+l] *= pivinv;
    for (ll=0; ll<n; ll++){
      if (ll != icol){
        idc = ll*n+icol;
        dum = Mat[idc];
        Mat[idc] = 0.;
        idc -= icol;
        for (l=0; l<n; l++) Mat[idc+l] -= Mat[idr+l]*dum;
      }
    }
  }
  for (l=n-1; l>=0; l--){
    int rl = indxr[l];
    int cl = indxc[l];
    if (rl != cl){
      for (k=0; k<n; k++){
        idr = k*n+rl;
        idc = k*n+cl;
        dum = Mat[idr];
        Mat[idr] = Mat[idc];
        Mat[idc] = dum;
      }
    }
  }

return;
}
