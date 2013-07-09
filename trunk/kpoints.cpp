#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "kpoints.h"

extern "C"{
#include "spglib.h"
}

#define MAXLINE 256
#define MAX(a,b) ((a) > (b) ? (a) : (b))

using namespace std;

/*------------------------------------------------------------------------------
 * Constructor is used as the driver to generate q-points
 *----------------------------------------------------------------------------*/
KPOINTS::KPOINTS(int type, CELL *cellin)
{
  memory = new Memory;
  cell = cellin;
  qr = NULL;

  if (type == 1){ // gamma point only
    nq   = 1;
    memory->create(q,nq, 3, "KPOINTS_KPOINTS:q");
    memory->create(w,nq,"KPOINTS_KPOINTS:w");
    q[0][0] = q[0][1] = q[0][2] = 0.;
    w[0] = 1.;

  } else if (type==3){ // line mode for dispersion
    get_line_q();

  } else {  // Monkhorst-Pack mesh
    int nx[3];
    char str[MAXLINE], *ptr;
    printf("\nPlease input the q-point mesh size [1 1 1]: ");
    if (cell->count_words(fgets(str,MAXLINE,stdin)) >= 3){
      ptr = strtok(str," \t\n\r\f");
      for (int i=0; i<3; i++){
        nx[i] = atoi(ptr);
        ptr = strtok(NULL," \t\n\r\f");
      }
    } else nx[0] = nx[1] = nx[2] = 1;

    for (int i=1; i<3; i++) if (nx[i] < 1) nx[i] = 1;
    for (int i=1; i<3; i++) if (cell->pbc[i] == 0) nx[i] = 1;
    printf("Your mesh size will be: %d x %d x %d\n", nx[0], nx[1], nx[2]);

    int npt = nx[0]*nx[1]*nx[2];
    get_ir_q(npt, nx, cell->natom, cell);
    printf("There are %d irreducible q-points in your mesh.\n", nq);

  }
return;
}

/*------------------------------------------------------------------------------
 * Free memories allocated
 *----------------------------------------------------------------------------*/
KPOINTS::~KPOINTS()
{
  memory->destroy(q);
  memory->destroy(w);
  memory->destroy(qr);

  cell = NULL;
  delete memory;
}

/*------------------------------------------------------------------------------
 * To get q-points for lines
 *----------------------------------------------------------------------------*/
void KPOINTS::get_line_q()
{
  // ask for q-points for lines
  int nline, *npt;
  double **qstr, **qend;
  char str[MAXLINE];

  nline = 0;
  while (nline < 1){
    nline = 1;
    printf("\nPlease input the total number of lines [1]: ");
    fgets(str,MAXLINE,stdin);
    char *ptr = strtok(str, " \t\n\r\f"); if (ptr) nline = atoi(ptr);
  }
  memory->create(qstr,nline, 3, "KPOINTS_interpolate:qstr");
  memory->create(qend,nline, 3, "KPOINTS_interpolate:qend");
  npt  = new int[nline];

  int npp = 10;
  for (int i=0; i<3; i++) qstr[0][i] = 0.;
  for (int i=0; i<nline; i++){
    printf("\nPlease input the start point of line %d [%g %g %g]: ", i+1, qstr[i][0], qstr[i][1],qstr[i][2]);
    if (cell->count_words(fgets(str,MAXLINE,stdin)) >= 3){
      qstr[i][0] = atof(strtok(str, " \t\n\r\f"));
      qstr[i][1] = atof(strtok(NULL," \t\n\r\f"));
      qstr[i][2] = atof(strtok(NULL," \t\n\r\f"));
    }

    do printf("Please input the final point of line %d: ", i+1);
    while(cell->count_words(fgets(str,MAXLINE,stdin)) < 3);
    qend[i][0] = atof(strtok(str, " \t\n\r\f"));
    qend[i][1] = atof(strtok(NULL," \t\n\r\f"));
    qend[i][2] = atof(strtok(NULL," \t\n\r\f"));

    printf("Please input the number of points in this line [%d]: ", npp);
    if (cell->count_words(fgets(str,MAXLINE,stdin)) > 0) npp = atoi(strtok(str, " \t\n\r\f"));
    npp = MAX(2,npp);
    npt[i] = npp;
    if (i < nline-1) for (int j=0; j<3; j++) qstr[i+1][j] = qend[i][j];
  }
  nq = 0;
  for (int i=0; i<nline; i++) nq += npt[i];
  memory->create(q, nq, 3, "KPOINTS_interpolate:q");
  memory->create(w, nq, "KPOINTS_interpolate:w");
  memory->create(qr,nq, "KPOINTS_interpolate:qr");
  
  // now to get the q-points
  double dq[3], r = 0., dr;
  nq = 0; // re-count, skip same points around the ends.
  for (int i=0; i<nline; i++){
    for (int k=0; k<3; k++) dq[k] = (qend[i][k] - qstr[i][k])/double(npt[i]-1);
    dr = sqrt(dq[0]*dq[0] + dq[1]*dq[1] + dq[2]*dq[2]);

    int istr = 0;
    if (i > 0){
      double dk[3];
      for (int k=0; k<3; k++) dk[k] = qstr[i][k]-qend[i-1][k];
      double diff = dk[0]*dk[0] + dk[1]*dk[1] + dk[2]*dk[2];
      if (diff < 1.e-8){istr = 1; r += dr;}
    }
    for (int j=istr; j<npt[i]; j++){
      for (int k=0; k<3; k++) q[nq][k] = qstr[i][k] + double(j)*dq[k];
      qr[nq] = r;
      r += dr;
      nq++;
    }
    r -= dr;
  }

  for (int i=0; i<nq; i++) w[i] = 1./double(nq);
  // free local memories
  memory->destroy(qstr);
  memory->destroy(qend);

return;
}

/*------------------------------------------------------------------------------
 * Private method to get the reciprocal vectors
 *----------------------------------------------------------------------------*/
void KPOINTS::writeq()
{
  // ask for file name
  char str[MAXLINE], *flag;
  while (1){
    printf("\nPlease input the file name to output q points:");
    fgets(str,MAXLINE,stdin);
    flag = strtok(str, " \t\n\r\f");
    if (flag) break;
  }

  int n = strlen(flag) + 1;
  char *fname = new char[n];
  strcpy(fname, flag);

  FILE *fp = fopen(fname, "w");
  fprintf(fp,"# index qx qy qz weight\n");
  for (int iq=0; iq<nq; iq++){
    fprintf(fp,"%d %lg %lg %lg %lg\n", iq+1, q[iq][0], q[iq][1], q[iq][2], w[iq]);
  }
  fclose(fp);
}

/*------------------------------------------------------------------------------
 * Private method to get the irreducible q-points and their weight
 *----------------------------------------------------------------------------*/
void KPOINTS::get_ir_q(const int npt, const int *nx, const int natom, CELL *cell)
{
  int mesh[3], shift[3], is_time_reversal = 0;
  int num_grid = npt, num_atom = natom;
  int grid_point[num_grid][3], map[num_grid], types[num_atom];
  double lattice[3][3], position[num_atom][3];
  
  for (int i=0; i<3; i++){ mesh[i] = nx[i]; shift[i] = 0; }
  for (int i=0; i<3; i++)
  for (int j=0; j<3; j++){
    lattice[i][j] = cell->axis[i][j];
  }
  for (int i=0; i<num_atom; i++){
    for (int j=0; j<3; j++) position[i][j] = cell->s[i][j];
    types[i] = cell->type[i];
  }
  double symprec = 1.e-5;

  nq = spg_get_ir_reciprocal_mesh(grid_point, map, num_grid,
                               mesh, shift, is_time_reversal,
                               lattice, position, types,
                               num_atom, symprec);

  memory->create(q,nq,3,"q");
  memory->create(w,nq,"w");
  int *iq2idx;
  memory->create(iq2idx, num_grid, "get_ir_q:iq2idx");
  int numq = 0;
  for (int i=0; i<num_grid; i++){
    int iq = map[i];
    if (iq == i) iq2idx[iq] = numq++;
  }
  for (int iq=0; iq<nq; iq++) w[iq] = 0.;
  numq = 0;
  for (int i=0; i<num_grid; i++){
    int iq = map[i];
    if (iq == i){ 
      q[numq][0] = double(grid_point[i][0])/double(mesh[0]);
      q[numq][1] = double(grid_point[i][1])/double(mesh[1]);
      q[numq][2] = double(grid_point[i][2])/double(mesh[2]);
      numq++;
    }
    w[iq2idx[iq]] += 1.;
  }

  memory->destroy(iq2idx);
  double wsum = 0.;
  for (int iq=0; iq<nq; iq++) wsum += w[iq];
  for (int iq=0; iq<nq; iq++) w[iq] /= wsum;

return;
}
