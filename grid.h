/* -*- c++ -*- ----------------------------------------------------------*/

#ifndef MPM_GRID_H
#define MPM_GRID_H

#include "pointers.h"
#include "grid.h"
#include "grid.h"
#include <vector>

class Grid : protected Pointers {
 public:
  double cellsize;                    // size of the square cells forming the grid
  vector< array<double, 3> > x;       // nodes' current position

  Grid(class MPM *);
  virtual ~Grid();
  void grow(int);
  void init(string);

protected:
  bigint nnodes; // number of particles
};

#endif
