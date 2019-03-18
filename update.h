/* -*- c++ -*- ----------------------------------------------------------*/

#ifndef MPM_UPDATE_H
#define MPM_UPDATE_H

#include "pointers.h"
#include <vector>

class Update : protected Pointers {
 public:
  double run_duration;            // GCG stop simulation if elapsed simulation time exceeds this.
  double elapsed_time_in_run;	  // elapsed simulation time for a single run;
  double dt;                      // timestep
  double dt_factor;               // timestep factor
  bigint ntimestep;               // current step
  int nsteps;                     // # of steps to run
  double atime;                   // simulation time at atime_step
  bigint atimestep;               // last timestep atime was updated
  bigint firststep,laststep;      // 1st & last step of this run
  bigint beginstep,endstep;       // 1st and last step of multiple runs
  int first_update;               // 0 before initial update, 1 after

  class Scheme *scheme;
  string scheme_style;

  class Method *method;
  string method_style;

  Update(class MPM *);
  ~Update();
  void set_dt_factor(vector<string>);
  void create_scheme(vector<string>);
  void new_scheme(vector<string>);
  void create_method(vector<string>);
  void new_method(vector<string>);
  void modify_method(vector<string>);
  void update_time();
protected:
  
};


#endif
