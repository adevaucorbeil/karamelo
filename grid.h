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
  int ncells;            // number of cells
  bigint nnodes;         // number of particles
  int nx;                // number of particles along x
  int ny;                // number of particles along y
  int nz;                // number of particles along z
  double cellsize;       // size of the square cells forming the grid

  // Eigen::Vector3d *C;            // connectivity matrix
  Eigen::Vector3d *x;            // nodes' current position
  Eigen::Vector3d *x0;           // nodes' reference position
  Eigen::Vector3d *v;            // nodes' velocity at time t
  Eigen::Vector3d *v_update;     // nodes' velocity at time t+dt
  Eigen::Vector3d *mb;           // nodes' external forces times the mass
  Eigen::Vector3d *f;            // nodes' internal forces

  // Eigen::Matrix3d *R;            // nodes' rotation matrix

  double *mass;              // nodes' current mass
  int *mask;                 // nodes' group mask
  int **ntype;               // node type in x, y, and z directions (False for an edge, True otherwise)

  Grid(class MPM *);
  virtual ~Grid();
  void grow(int);
  void setup(string);
  void init(double*, double*);

  void update_grid_velocities();
  void update_grid_positions();
};

#endif
