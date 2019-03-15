/* -*- c++ -*- ----------------------------------------------------------*/

#ifndef MPM_FIX_H
#define MPM_FIX_H

#include "pointers.h"
#include <vector>

class Fix : protected Pointers {
 public:
  string id;
  int igroup, groupbit;
  int mask;

  Fix(class MPM *, vector<string>);
  virtual ~Fix() {};
  virtual void init() = 0;
  virtual void setup() = 0;
  virtual void setmask() = 0;

  
  virtual void initial_integrate() = 0;
  virtual void post_update_grid_state() = 0;
  virtual void post_grid_to_point() = 0;
  virtual void post_advance_particles() = 0;
  virtual void post_velocities_to_grid() = 0;
  virtual void final_integrate() = 0;
};

namespace FixConst {
  static const int INITIAL_INTEGRATE =       1<<0;
  static const int POST_UPDATE_GRID_STATE =  1<<1;
  static const int POST_GRID_TO_POINT =      1<<2;
  static const int POST_ADVANCE_PARTICLES =  1<<3;
  static const int POST_VELOCITIES_TO_GRID = 1<<4;
  static const int FINAL_INTEGRATE =         1<<5;
}
#endif
