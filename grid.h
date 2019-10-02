/* -*- c++ -*- ----------------------------------------------------------*/

#ifndef MPM_GRID_H
#define MPM_GRID_H

#include "pointers.h"
#include "grid.h"
#include "grid.h"
#include <vector>
#include <Eigen/Eigen>
#include <map>

using namespace Eigen;

struct Point {
  tagint tag;
  double x[3];
  int owner;
};

class Grid : protected Pointers {
 public:
  int ncells;            // number of cells
  bigint nnodes;         // total number of nodes in the domain
  bigint nnodes_local;   // number of nodes (in this CPU)
  bigint nnodes_ghost;   // number of ghost nodes (in this CPU)
  tagint *ntag;          // unique identifier for nodes in the system.
  map<int, int> map_ntag;// map_ntag[ntag[i]] = i;

  int nx;                // number of nodes along x on this CPU
  int ny;                // number of nodes along y on this CPU
  int nz;                // number of nodes along z on this CPU

  int nx_global;         // number of nodes along x on all CPUs
  int ny_global;         // number of nodes along y on all CPUs
  int nz_global;         // number of nodes along z on all CPUs

  int nshared;           // number of nodes that are shared (ghosts in other CPUs
  vector<int> shared;    // position of all shared nodes

  int *nowner;          // which CPU owns each node (universe->me for local nodes, other CPU for ghost nodes

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

  MPI_Datatype Pointtype;    // MPI type for struct Point

  Grid(class MPM *);
  virtual ~Grid();
  void grow(int);
  void grow_ghosts(int);
  void setup(string);
  void init(double*, double*);

  void reduce_ghost_nodes(bool only_v = false);
  void update_grid_velocities();
  void update_grid_positions();
};

#endif
