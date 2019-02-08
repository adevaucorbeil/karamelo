/* -*- c++ -*- ----------------------------------------------------------*/

#ifndef MPM_GRID_H
#define MPM_GRID_H

#include "pointers.h"
#include "grid.h"
#include "grid.h"
#include <vector>

class Grid : protected Pointers {
 public:
  double cellsize;       // size of the square cells forming the grid
  double **x;            // nodes' current position
  double **v;            // nodes' velocity at time t
  double **v_update;     // nodes' velocity at time t+dt
  double **b;            // nodes' external forces
  double **f;            // nodes' internal forces

  Grid(class MPM *);
  virtual ~Grid();
  void grow(int);
  void init(string);

protected:
  bigint nnodes; // number of particles
};

#endif
