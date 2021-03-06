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

FixStyle(contact/minimize_penetration, FixContactMinPenetration)

#else

#ifndef MPM_FIX_CONTACT_MIN_PENETRATION_H
#define MPM_FIX_CONTACT_MIN_PENETRATION_H

#include "fix.h"
#include "var.h"
#include <vector>

class FixContactMinPenetration : public Fix {
public:
  FixContactMinPenetration(class MPM *, vector<string>);
  ~FixContactMinPenetration();
  void setmask();
  void init();
  void setup();

  void initial_integrate();
  void post_particles_to_grid(){};
  void post_update_grid_state(){};
  void post_grid_to_point(){};
  void post_advance_particles() {};
  void post_velocities_to_grid(){};
  void final_integrate(){};

  void write_restart(ofstream *);
  void read_restart(ifstream *);

private:
  string usage = "Usage: fix(fix-ID, contact/minimize_penetration, solid1, solid2, mu)\n";
  int Nargs = 5;
  int solid1, solid2;
  double mu;    // Friction coefficient
};

#endif
#endif

