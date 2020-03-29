/* -*- c++ -*- ----------------------------------------------------------*/

#ifdef FIX_CLASS

FixStyle(velocity_particles,FixVelocityParticles)

#else

#ifndef MPM_FIX_VELOCITY_PARTICLES_H
#define MPM_FIX_VELOCITY_PARTICLES_H

#include "fix.h"
#include "var.h"
#include <Eigen/Eigen>
#include <vector>

class FixVelocityParticles : public Fix {
 public:
  FixVelocityParticles(class MPM *, vector<string>);
  ~FixVelocityParticles();
  void setmask();
  void init();
  void setup();
  
  void initial_integrate();
  void post_particles_to_grid() {};
  void post_update_grid_state() {};
  void post_grid_to_point() {};
  void post_advance_particles();
  void post_velocities_to_grid() {};
  void final_integrate() {};

private:
  //class Var xvalue, yvalue, zvalue;    // Set velocities in x, y, and z directions.
  bool xset, yset, zset;                 // Does the fix set the x, y, and z velocities of the group?
  int xpos, ypos, zpos;                  // Positions of x, y and z in the argument list (args)
  vector<string> args_previous_step;

  vector<Eigen::Vector3d> xold;          // particles' old position
};

#endif
#endif

