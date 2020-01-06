/* -*- c++ -*- ----------------------------------------------------------
 *
 *                    ***       Karamelo       ***
 *               Parallel Material Point Method Simulator
 * 
 * Copyright (2019) Alban de Vaucorbeil, alban.devaucorbeil@monash.edu
 * Materials Science and Engineering, Monash University
 * Clayton VIC 3800, Australia

 * This software is distributed under the GNU General Public License.
 *
 * ----------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(body_force,FixBodyforce)

#else

#ifndef MPM_FIX_BODY_FORCE_H
#define MPM_FIX_BODY_FORCE_H

#include "fix.h"
#include "var.h"
#include <vector>

class FixBodyforce : public Fix {
 public:
  FixBodyforce(class MPM *, vector<string>);
  ~FixBodyforce();
  void setmask();
  void init();
  void setup();
  
  void initial_integrate() {};
  void post_particles_to_grid();
  void post_update_grid_state() {};
  void post_grid_to_point() {};
  void post_advance_particles() {};
  void post_velocities_to_grid() {};
  void final_integrate() {};

private:
  class Var xvalue, yvalue, zvalue;    // Set force in x, y, and z directions.
  bool xset, yset, zset;               // Does the fix set the x, y, and z forces of the group?
};

#endif
#endif

