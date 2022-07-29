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

#ifdef DUMP_CLASS

DumpStyle(grid/bin,DumpGridBin)

#else

#ifndef MPM_DUMP_GRID_BIN_H
#define MPM_DUMP_GRID_BIN_H

#include <deque>
#include <dump.h>
#include <grid.h>
#include <thread>
#include <utility>

class DumpGridBin : public Dump {
  Kokkos::View<tagint*>::HostMirror ntag;   ///< unique identifier for nodes in the system.
  
  Kokkos::View<Vector3d*>::HostMirror x;            ///< nodes' current position
  Kokkos::View<Vector3d*>::HostMirror v;            ///< nodes' velocity at time t
  Kokkos::View<Vector3d*>::HostMirror mb;           ///< nodes' external forces times the mass

  Kokkos::View<float*>::HostMirror mass;              ///< nodes' current mass
  Kokkos::View<int*>::HostMirror mask;                 ///< nodes' group mask
  Kokkos::View<bool*>::HostMirror rigid;               ///< are the nodes in the area of influence of a rigid body?
  Kokkos::View<Vector3i*>::HostMirror ntype;      ///< node type in x, y, and z directions (False for an edge, True otherwise)

  Kokkos::View<float*>::HostMirror T;                 ///< nodes' temperature at time t

 public:
  DumpGridBin(MPM *, vector<string>);
  ~DumpGridBin();

  void write();
 protected:
  vector<string> known_var = {"x", "y", "z",
			      "vx", "vy", "vz",
			      "bx", "by", "bz",
			      "mass", "mask",
			      "rigid", "T",
			      "ntypex", "ntypey", "ntypez"};
private:
  const string MAGIC_STRING = "DUMPCUSTOM";
  const int FORMAT_REVISION = 0x0002;
  const int ENDIAN = 0x0001;

  deque<pair<thread, vector<double>>>  threads;        ///< Pair storing the threads and the buffer

  void write_to_file(bigint, string, bigint, bigint);
};

#endif
#endif
