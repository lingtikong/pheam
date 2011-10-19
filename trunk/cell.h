#ifndef CELL_H
#define CELL_H

#include "memory.h"
#include <math.h>

class CELL {
public:
  CELL(void);       // to read in cell information
  ~CELL();          // to free memory

  Memory *memory;

  int natom, ntype; // # of atoms and atom types
  int pbc[3];       // 0, free; 1, periodic boundary condition
  double **axis, **invaxis; // cell vector matrix

  int *type;        // atomic type
  double **x, **s;  // atom positions
  char **elements;  // element names
  char *title;      // title of the cell

  int index(const char *);// index of elements in type
  double veclen(double *); // the vector length; vector should be in fractional
  double veclen2(double *);// square of the vector length; vector should be in fractional
  double VecLen(double *); // the vector length; vector should be in cartesian
  double VecLen2(double *);// square of the vector length; vector should be in cartesian
  double VecAng(double *, double *); // the cosine between two vectors, should be in cartesian

  void display();   // display the cell info
  void GaussJordan(const int, const double *, double *); // To get the inverse of a square matrix

  void dir2car(double *, double *); // convert a fractional vector into a cartesian one

  int count_words(const char *);    // method to count the number of words in a string

private:
  void dir2car();   // to convert atomic coordinate from directional to cartesian
  void car2dir();   // to convert atomic coordinate from cartesian to directional
};
#endif
