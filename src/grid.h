/* -*- c++ -*- ----------------------------------------------------------
 *
 *                    ***       Karamelo       ***
 *               Parallel Material Point Method Simulator
 * 
 * Copyright (2019) Alban de Vaucorbeil, alban.devaucorbeil@monash.edu
 * Materials Science and Engineering, Monash University
 * Clayton VIC 3800, Australia

 * This software is distributed under the GNU General Public License.
 *
 * ----------------------------------------------------------------------- */

#ifndef MPM_GRID_H
#define MPM_GRID_H

#include "pointers.h"
#include "grid.h"
#include "grid.h"
#include <vector>
#include <Eigen/Eigen>
#include <map>
#include <array>

using namespace Eigen;

struct Point {
  int owner;
  tagint tag;
  double x[3];
  int ntype[3];
};

class Grid : protected Pointers {
 public:
  int ncells;            // number of cells
  bigint nnodes;         // total number of nodes in the domain
  bigint nnodes_local;   // number of nodes (in this CPU)
  bigint nnodes_ghost;   // number of ghost nodes (in this CPU)
  vector<tagint> ntag;   // unique identifier for nodes in the system.
  map<int, int> map_ntag;// map_ntag[ntag[i]] = i;

  int nx;                // number of nodes along x on this CPU
  int ny;                // number of nodes along y on this CPU
  int nz;                // number of nodes along z on this CPU

  int nx_global;         // number of nodes along x on all CPUs
  int ny_global;         // number of nodes along y on all CPUs
  int nz_global;         // number of nodes along z on all CPUs

  int nshared;           // number of nodes that are shared (ghosts in other CPUs
  vector<tagint> shared; // tag of all shared nodes
  map<int, vector<tagint>> dest_nshared;   // for each CPU, list the tags of shared nodes
  map<int, vector<tagint>> origin_nshared; // for each CPU, list the tags of ghost nodes

  vector<int> nowner;    // which CPU owns each node (universe->me for local nodes, other CPU for ghost nodes

  double cellsize;       // size of the square cells forming the grid

  vector<Eigen::Vector3d> x;            // nodes' current position
  vector<Eigen::Vector3d> x0;           // nodes' reference position
  vector<Eigen::Vector3d> v;            // nodes' velocity at time t
  vector<Eigen::Vector3d> v_update;     // nodes' velocity at time t+dt
  vector<Eigen::Vector3d> mb;           // nodes' external forces times the mass
  vector<Eigen::Vector3d> f;            // nodes' internal forces

  vector<double> mass;              // nodes' current mass
  vector<int> mask;                 // nodes' group mask
  vector<bool> rigid;               // are the nodes in the area of influence of a rigid body?
  vector<array<int, 3>> ntype;      // node type in x, y, and z directions (False for an edge, True otherwise)

  MPI_Datatype Pointtype;    // MPI type for struct Point

  Grid(class MPM *);
  virtual ~Grid();
  void grow(int);
  void grow_ghosts(int);
  void setup(string);
  void init(double*, double*);

  void reduce_mass_ghost_nodes();
  void reduce_mass_ghost_nodes_old();
  void reduce_ghost_nodes(bool only_v = false);
  void reduce_ghost_nodes_old(bool only_v = false);
  void update_grid_velocities();
  void update_grid_positions();
};

#endif