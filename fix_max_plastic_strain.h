/* -*- c++ -*- ----------------------------------------------------------*/

#ifdef FIX_CLASS

FixStyle(max_plastic_strain, FixMaxPlasticStrain)

#else

#ifndef MPM_FIX_MAX_PLASTIC_STRAIN_H
#define MPM_FIX_MAX_PLASTIC_STRAIN_H

#include "fix.h"
#include "var.h"
#include <vector>

class FixMaxPlasticStrain : public Fix
{
public:
  FixMaxPlasticStrain(class MPM *, vector<string>);
  ~FixMaxPlasticStrain();
  void setmask();
  void init();
  void setup();

  void initial_integrate(){};
  void post_particles_to_grid(){};
  void post_update_grid_state(){};
  void post_grid_to_point(){};
  void post_advance_particles(){};
  void post_velocities_to_grid(){};
  void final_integrate();

private:
  // class Var xvalue, yvalue, zvalue;    // Set force in x, y, and z
  // directions. bool xset, yset, zset;               // Does the fix set the x,
  // y, and z forces of the group?
};

#endif
#endif
