/* -*- c++ -*- ----------------------------------------------------------*/

#ifdef FIX_CLASS

FixStyle(kinetic_energy,FixKineticEnergy)

#else

#ifndef MPM_FIX_KINETIC_ENERGY_H
#define MPM_FIX_KINETIC_ENERGY_H

#include "fix.h"
#include "var.h"
#include <vector>

class FixKineticEnergy : public Fix {
 public:
  FixKineticEnergy(class MPM *, vector<string>);
  ~FixKineticEnergy();
  void setmask();
  void init();
  void setup();
  
  void initial_integrate() {};
  void post_particles_to_grid() {};
  void post_update_grid_state() {};
  void post_grid_to_point() {};
  void post_advance_particles() {};
  void post_velocities_to_grid() {};
  void final_integrate();

private:
  // class Var xvalue, yvalue, zvalue;    // Set force in x, y, and z directions.
  // bool xset, yset, zset;               // Does the fix set the x, y, and z forces of the group?
};

#endif
#endif

