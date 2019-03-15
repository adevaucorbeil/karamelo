/* -*- c++ -*- ----------------------------------------------------------*/

#ifdef FIX_CLASS

FixStyle(velocity,FixVelocity)

#else

#ifndef MPM_FIX_VELOCITY_H
#define MPM_FIX_VELOCITY_H

#include "fix.h"
#include "var.h"
#include <vector>

class FixVelocity : public Fix {
 public:
  FixVelocity(class MPM *, vector<string>);
  ~FixVelocity();
  // int setmask();
  void init();
  void setup();
  // void min_setup(int);
  // void initial_integrate(int);
  // void post_force(int);
  // double compute_vector(int);
  // double memory_usage(); 

private:
  class Var xvalue, yvalue, zvalue;    // Set velocities in x, y, and z directions.
  bool xset, yset, zset;               // Does the fix set the x, y, and z velocities of the group?
};

#endif
#endif

