/* -*- c++ -*- ----------------------------------------------------------*/

#ifdef FIX_CLASS

FixStyle(initial_velocity_particles,FixInitialVelocityParticles)

#else

#ifndef MPM_FIX_INITIAL_VELOCITY_PARTICLES_H
#define MPM_FIX_INITIAL_VELOCITY_PARTICLES_H

#include "fix.h"
#include "var.h"
#include <vector>

class FixInitialVelocityParticles : public Fix {
 public:
  FixInitialVelocityParticles(class MPM *, vector<string>);
  ~FixInitialVelocityParticles();
  void setmask();
  void init();
  void setup();
  
  void initial_integrate();
  void post_particles_to_grid() {};
  void post_update_grid_state() {};
  void post_grid_to_point() {};
  void post_advance_particles() {};
  void post_velocities_to_grid() {};
  void final_integrate() {};

private:
  class Var xvalue, yvalue, zvalue;    // Set velocities in x, y, and z directions.
  bool xset, yset, zset;               // Does the fix set the x, y, and z velocities of the group?
};

#endif
#endif

