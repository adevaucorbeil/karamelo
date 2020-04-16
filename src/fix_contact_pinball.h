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

FixStyle(contact/pinball, FixContactPinball)

#else

#ifndef MPM_FIX_CONTACT_PINBALL_H
#define MPM_FIX_CONTACT_PINBALL_H

#include "fix.h"
#include "var.h"
#include <vector>

class FixContactPinball : public Fix {
public:
  FixContactPinball(class MPM *, vector<string>);
  ~FixContactPinball();
  void setmask();
  void init();
  void setup();

  void initial_integrate();
  void post_particles_to_grid(){};
  void post_update_grid_state(){};
  void post_grid_to_point(){};
  void post_advance_particles(){};
  void post_velocities_to_grid(){};
  void final_integrate(){};

private:
  string usage = "Usage: fix(fix-ID, contact/pinball, solid1, solid2, K)\n";
  int Nargs = 5;
  int solid1, solid2;
  double K;
};

#endif
#endif

