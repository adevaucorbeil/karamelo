/* -*- c++ -*- ----------------------------------------------------------*/


#ifndef LMP_UNIVERSE_H
#define LMP_UNIVERSE_H

#include "pointers.h"

using namespace std;

class Universe : protected Pointers {
 public:
  const char *version;    // Karamelo version string = date

  MPI_Comm uworld;        // communicator for entire universe
  int me,nprocs;          // my place in universe

  Universe(class MPM *, MPI_Comm);
  ~Universe();

  int procgrid[3];                  // procs assigned in each dim of 3d grid
  int myloc[3];                     // which proc I am in each dim
  int procneigh[3][2];              // my 6 neighboring procs, 0/1 = left/right

  void set_proc_grid();             // setup 3d grid of procs
};

#endif
