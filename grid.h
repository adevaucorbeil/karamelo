/* -*- c++ -*- ----------------------------------------------------------*/

#ifndef MPM_GRID_H
#define MPM_GRID_H

#include "pointers.h"
#include "grid.h"
#include "grid.h"
#include <vector>
#include <Eigen/Eigen>

using namespace Eigen;

class Grid : protected Pointers {
 public:
  double cellsize;       // size of the square cells forming the grid

  Eigen::Vector3d *x;            // nodes' current position
  Eigen::Vector3d *v;            // nodes' velocity at time t
  Eigen::Vector3d *v_update;     // nodes' velocity at time t+dt
  Eigen::Vector3d *b;            // nodes' external forces
  Eigen::Vector3d *f;            // nodes' internal forces


  Grid(class MPM *);
  virtual ~Grid();
  void grow(int);
  void init(string);

protected:
  bigint nnodes; // number of particles
};

#endif
