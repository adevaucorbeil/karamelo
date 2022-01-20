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

#include <pointers.h>
#include <grid.h>
#include <grid.h>
#include <vector>
#include <matrix.h>
#include <unordered_map>
#include <map>
#include <array>



/*! This structure is used to duplicate a grid point to another CPU.
 *  
 * Each CPU create the series of grid points that lie in their respective domains.
 * To insure the continuity of the simulation domain, each CPU need to be aware of
 * of the grid points in the vicinity of its domain. 
 * These points are communicated by the CPU whose domain they lie onto (owner).
 * Each point has a unique identification number (tag), and a position vector (x).
 * Moreover, it has a type that is used when B-splines and Bernstein quadratic shape 
 * functions are used. This type changes in function of the position of the point with
 * respect to the boundaries of the whole domain, or the quadratic cell.
 */
struct Point {
  int owner;    ///< ID of the CPU who created the point
  tagint tag;   ///< Unique identification number of the point
  double x[3];  ///< Position vector of the point
  int ntype[3]; ///< Type of the point (0, 1, or 2)
};


/*! This class holds all the variables and independent functions related the background grid(s).
 * 
 * This class does not have as members functions related to the Particle to Grid (P2G) and Grid to Particle (G2P) steps of the MPM algorithm.
 * This is due to the necessity to access the mass, velocity, body forces and stresses of particles held in the Solid class.
 * Instead, these functions are members of the Solid Class since this class stores a pointer to the grid it uses. 
 * Either a local grid, as used in the Total Lagrangian MPM, or the global grid, as used in the Updated Lagrangian MPM.
 * This class has Pointers as parent in order to be accessible from anywhere within the MPM class.
 */
class Grid : protected Pointers {
 public:
  int ncells;            ///< number of cells
  bigint nnodes;         ///< total number of nodes in the domain
  bigint nnodes_local;   ///< number of nodes (in this CPU)
  bigint nnodes_ghost;   ///< number of ghost nodes (in this CPU)
  vector<tagint> ntag;   ///< unique identifier for nodes in the system.
  vector<tagint> map_ntag;  ///< map_ntag[ntag[i]] = i;

  int nx;                ///< number of nodes along x on this CPU
  int ny;                ///< number of nodes along y on this CPU
  int nz;                ///< number of nodes along z on this CPU

  int nx_global;         ///< number of nodes along x on all CPUs
  int ny_global;         ///< number of nodes along y on all CPUs
  int nz_global;         ///< number of nodes along z on all CPUs

  int nshared;           ///< number of nodes that are shared (ghosts in other CPUs)
  vector<tagint> shared; ///< tag of all shared nodes
  map<int, vector<tagint>> dest_nshared;   ///< for each CPU, list the tags of shared nodes
  map<int, vector<tagint>> origin_nshared; ///< for each CPU, list the tags of ghost nodes

  vector<int> nowner;    ///< which CPU owns each node (universe->me for local nodes, other CPU for ghost nodes

  double cellsize;       ///< size of the square cells forming the grid

  vector<Vector3d> x;            ///< nodes' current position
  vector<Vector3d> x0;           ///< nodes' position in the reference coordinate system
  vector<Vector3d> v;            ///< nodes' velocity at time t
  vector<Vector3d> v_update;     ///< nodes' velocity at time t+dt
  vector<Vector3d> mb;           ///< nodes' external forces times the mass
  vector<Vector3d> f;            ///< nodes' internal forces

  vector<double> mass;              ///< nodes' current mass
  vector<int> mask;                 ///< nodes' group mask
  vector<bool> rigid;               ///< are the nodes in the area of influence of a rigid body?
  vector<array<int, 3>> ntype;      ///< node type in x, y, and z directions (False for an edge, True otherwise)

  vector<double> T;                 ///< nodes' temperature at time t
  vector<double> T_update;          ///< nodes' temperature at time t+dt
  vector<double> Qext;              ///< nodes' external thermal driving force
  vector<double> Qint;              ///< nodes' internal thermal driving force

  MPI_Datatype Pointtype;           ///< MPI type for struct Point

  Grid(class MPM *);
  virtual ~Grid();
  void grow(int);              ///< Allocate memory for the vectors used for local nodes or resize them  
  void grow_ghosts(int);       ///< Allocate memory for the vectors used for ghost nodes or resize them
  void setup(string);
  void init(double*, double*); ///< Create the array of nodes. Give them their position, tag, and type

  void reduce_mass_ghost_nodes();                  ///< Reduce the mass of all the ghost nodes from that computed on each CPU.
  void reduce_mass_ghost_nodes_old();              ///< Deprecated
  void reduce_rigid_ghost_nodes();                 ///< Reduce the rigid bool of all the ghost nodes from that computed on each CPU.
  void reduce_ghost_nodes(bool reduce_v, bool reduce_forces, bool temp = false);    ///< Reduce the force and velocities of all the ghost nodes from that computed on each CPU.
  void reduce_ghost_nodes_old(bool only_v = false, bool temp = false);    ///< Deprecated
  void update_grid_velocities();                   ///< Determine the temporary grid velocities \f$\tilde{v}_{n}\f$. 
  void update_grid_positions();                    ///< Determine the new position of the grid nodes.
  void update_grid_temperature();                  ///< Determine the temporary grid temperature \f$\tilde{T}_{n}\f$.

  void reset_mass();
  void reset_velocity();
  void reset_forces();
};

#endif
