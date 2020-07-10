/* -*- c++ -*- ----------------------------------------------------------*/

#ifdef FIX_CLASS

FixStyle(velocity_particles_exact_integration,FixVelocityParticlesExactIntegration)

#else

#ifndef MPM_FIX_VELOCITY_PARTICLES_EXACT_INTEGRATION_H
#define MPM_FIX_VELOCITY_PARTICLES_EXACT_INTEGRATION_H

#include "fix.h"
#include "var.h"
#include <Eigen/Eigen>
#include <vector>

class FixVelocityParticlesExactIntegration : public Fix {
 public:
  FixVelocityParticlesExactIntegration(class MPM *, vector<string>);
  ~FixVelocityParticlesExactIntegration();
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
  int xpos = 3;
  int vxpos = 4;
  int ypos = 5;
  int vypos = 6;
  int zpos = 7;
  int vzpos = 8;
  vector<string> args_previous_step, args_next_step;

  vector<Eigen::Vector3d> xold;          // particles' old position

  map<int, string> usage = {
      {1, "Usage: fix(fix-ID, velocity_particles_exact_integration, group-ID, "
          "x, vx)\n"},
      {2, "Usage: fix(fix-ID, velocity_particles_exact_integration, group-ID, "
          "x, vx, y, "
          "vy)\n"},
      {3, "Usage: fix(fix-ID, velocity_particles_exact_integration, group-ID, "
          "x, vx, y, "
          "vy, z, vz)\n"}};
  map<int, int> Nargs = { {1,5}, {2,7}, {3,9}} ;
};

#endif
#endif

