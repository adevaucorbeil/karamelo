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
#include <vector>
#include <matrix.h>
#include <unordered_map>
#include <map>
#include <array>
#include <Kokkos_Core.hpp>

#ifdef KOKKOS_ENABLE_CUDA
using MemorySpace = Kokkos::CudaUVMSpace;
#else
using MemorySpace = Kokkos::HostSpace;
#endif

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
class Grid : public Pointers {
 public:
  int ncells;            ///< number of cells
  bigint nnodes;         ///< total number of nodes in the domain
  bigint nnodes_local;   ///< number of nodes (in this CPU)
  bigint nnodes_ghost;   ///< number of ghost nodes (in this CPU)
  Kokkos::View<tagint*, MemorySpace> ntag;   ///< unique identifier for nodes in the system.
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

  Kokkos::View<int*, MemorySpace> nowner;    ///< which CPU owns each node (universe->me for local nodes, other CPU for ghost nodes

  double cellsize;       ///< size of the square cells forming the grid

  Kokkos::View<Vector3d*, MemorySpace> x;            ///< nodes' current position
  Kokkos::View<Vector3d*, MemorySpace> x0;           ///< nodes' position in the reference coordinate system
  Kokkos::View<Vector3d*, MemorySpace> v;            ///< nodes' velocity at time t
  Kokkos::View<Vector3d*, MemorySpace> v_update;     ///< nodes' velocity at time t+dt
  Kokkos::View<Vector3d*, MemorySpace> mb;           ///< nodes' external forces times the mass
  Kokkos::View<Vector3d*, MemorySpace> f;            ///< nodes' internal forces

  Kokkos::View<double*, MemorySpace> mass;              ///< nodes' current mass
  Kokkos::View<int*, MemorySpace> mask;                 ///< nodes' group mask
  Kokkos::View<bool*, MemorySpace> rigid;               ///< are the nodes in the area of influence of a rigid body?
  Kokkos::View<Vector3i*, MemorySpace> ntype;      ///< node type in x, y, and z directions (False for an edge, True otherwise)

  Kokkos::View<double*, MemorySpace> T;                 ///< nodes' temperature at time t
  Kokkos::View<double*, MemorySpace> T_update;          ///< nodes' temperature at time t+dt
  Kokkos::View<double*, MemorySpace> Qext;              ///< nodes' external thermal driving force
  Kokkos::View<double*, MemorySpace> Qint;              ///< nodes' internal thermal driving force

  Grid(class MPM *);
  void grow(int);              ///< Allocate memory for the vectors used for local nodes or resize them  
  void setup(string);
  void init(double*, double*); ///< Create the array of nodes. Give them their position, tag, and type

  void reduce_mass_ghost_nodes();                  ///< Reduce the mass of all the ghost nodes from that computed on each CPU.
  void reduce_mass_ghost_nodes_old();              ///< Deprecated
  void reduce_rigid_ghost_nodes();                 ///< Reduce the rigid bool of all the ghost nodes from that computed on each CPU.
  void reduce_ghost_nodes(bool reduce_v, bool reduce_forces, bool temp = false);    ///< Reduce the force and velocities of all the ghost nodes from that computed on each CPU.
  void reduce_ghost_nodes_old(bool only_v = false, bool temp = false);    ///< Deprecated
};

#endif
