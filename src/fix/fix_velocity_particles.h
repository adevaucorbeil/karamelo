/* -*- c++ -*- ----------------------------------------------------------*/

#ifdef FIX_CLASS

FixStyle(velocity_particles,FixVelocityParticles)

#else

#ifndef MPM_FIX_VELOCITY_PARTICLES_H
#define MPM_FIX_VELOCITY_PARTICLES_H

#include <fix.h>
#include <var.h>
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

  void write_restart(ofstream *);
  void read_restart(ifstream *);

private:
  const map<int, string> usage = {
      {1, "Usage: fix(fix-ID, velocity_particles, group, vx)\n"},
      {2, "Usage: fix(fix-ID, velocity_particles, group, vx, vy)\n"},
      {3, "Usage: fix(fix-ID, velocity_particles, group, vx, vy, vz)\n"}};
  const map<int, int> Nargs = {{1, 4}, {2, 5}, {3, 6}};

  class Var xvalue, yvalue, zvalue;                  //< Velocities in x, y, and z directions.
  class Var xprevvalue, yprevvalue, zprevvalue;      //< Velocities in x, y, and z directions from previous time step.
  bool xset, yset, zset;                             //< Does the fix set the x, y, and z velocities of the group?

  vector<Eigen::Vector3d> xold;          // particles' old position
};

#endif
#endif

