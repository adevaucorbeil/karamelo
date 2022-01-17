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


#ifndef LMP_UNIVERSE_H
#define LMP_UNIVERSE_H

#include "pointers.h"
#include <vector>

using namespace std;

/*! Sets up partitions of processors so that multiple simulations can be run, 
 * each on a subset of the processors allocated for a run, e.g. by the mpirun command
 */
class Universe : protected Pointers {
 public:
  MPI_Comm uworld;        ///< MPI communicator for the entire universe
  int me,nprocs;          ///< My place (as a proc) in the universe

  Universe(class MPM *, MPI_Comm);
  ~Universe();

  int procgrid[3];                  ///< procs assigned in each dim of 3d grid
  int myloc[3];                     ///< which proc I am in each dim
  int procneigh[3][2];              ///< my 6 neighboring procs, 0/1 = left/right

  vector<array<int, 3>> sendnrecv;  ///< stores the send and receive pattern to follow. The vector's length corresponds to the number of steps. For a given step, there is an array of two ints, the first int is either 0 (receive) or 1 (send), the second it the rank of the CPU to receive or send to.

  void set_proc_grid();             ///< setup 3d grid of procs
};

#endif
