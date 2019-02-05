/* -*- c++ -*- ----------------------------------------------------------*/

#ifndef MPM_GRID_H
#define MPM_GRID_H

#include "pointers.h"
#include "grid.h"
#include <vector>

class Grid : protected Pointers {
 public:
  vector< array<double, 3> > x;       // nodes' current position

  Grid(class MPM *);
  virtual ~Grid();
  void grow(int);

protected:
  bigint nnodes; // number of particles
};

#endif
