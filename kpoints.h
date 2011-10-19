#ifndef KPOINTS_H
#define KPOINTS_H

#include "memory.h"
#include "cell.h"

class KPOINTS {
public:
  KPOINTS(int, CELL *);
  ~KPOINTS();

  int nq;
  double **q, *w, *qr;  // q points and their weight

  void writeq();

private:
  Memory *memory;
  CELL *cell;

  void get_ir_q(const int, const int *, const int, CELL *);
  void get_line_q();
};
#endif
