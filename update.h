/* -*- c++ -*- ----------------------------------------------------------*/

#ifndef MPM_UPDATE_H
#define MPM_UPDATE_H

#include "pointers.h"

class Update : protected Pointers {
 public:
  double run_duration;            // GCG stop simulation if elapsed simulation time exceeds this.
  double elapsed_time_in_run;	  // elapsed simulation time for a single run;
  double dt;                      // timestep
  bigint ntimestep;               // current step
  int nsteps;                     // # of steps to run
  double atime;                   // simulation time at atime_step
  bigint atimestep;               // last timestep atime was updated
  bigint firststep,laststep;      // 1st & last step of this run
  bigint beginstep,endstep;       // 1st and last step of multiple runs
  int first_update;               // 0 before initial update, 1 after

  Update(class MPM *);
  ~Update();

};


#endif
