#ifndef EAM_H
#define EAM_H

#include "memory.h"

class EAM {
public:
  
  EAM();
  ~EAM();

  int initialized;  // 0, not initialized; 1, initialized.
  int ntype;        // # of atomic types in the EAM file
  char **elements;  // element names of the EAM
  double *mass;     // masses of the elements
  double dr, drho, rcut, rcutsq; // increment for distance and electronic density, as well as cutoffs

  Memory *memory;   // class memory to create and destroy arrays.

  virtual double F(double, int){ return 0.; }          // the embedded function
  virtual double Rho(double, int, int){ return 0.; }   // the electronic density function
  virtual double Phi(double, int, int){ return 0.; }   // the pair interaction function

  virtual double Fp(double, int){ return 0.; }         // F'(rho)
  virtual double Rhop(double, int, int){ return 0.; }  // f'(r)
  virtual double Phip(double, int, int){ return 0.; }  // P'(r)

  virtual double Fpp(double, int){ return 0.; }        // F"(rho)
  virtual double Rhopp(double, int, int){ return 0.; } // f"(r)
  virtual double Phipp(double, int, int){ return 0.; } // P"(r)

 // double LagrangeInterpolate(double, double *); // method to do four point Lagrange interpolation
  int index(const char *);        // method to find the index of an element type in EAM
  int count_words(const char *);  // method to count the number of words in a string

  void grab(FILE *, int, double *);

  void check_eam();

private:
  void lateng(double, double *, double *, int, int);
};
#endif
