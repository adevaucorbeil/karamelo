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

FixStyle(initial_stress,FixInitialStress)

#else

#ifndef MPM_FIX_INITIAL_STRESS_H
#define MPM_FIX_INITIAL_STRESS_H

#include <fix.h>
#include <var.h>
#include <vector>

class FixInitialStress : public Fix {
 public:
  FixInitialStress(class MPM *, vector<string>);
  ~FixInitialStress();
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

  void write_restart(ofstream *) {};
  void read_restart(ifstream *) {};

private:
  string usage = "Usage: fix(fix-ID, initial_stress, group-ID, sigma_xx, sigma_yy, sigma_zz, sigma_yz, sigma_xz, sigma_xy)\n";
  int Nargs = 9;

  class Var s_value[6];
  bool s_set[6];
};

#endif
#endif

