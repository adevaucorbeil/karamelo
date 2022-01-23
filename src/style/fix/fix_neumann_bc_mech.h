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

FixStyle(neumann_bc_mech,FixNeumannBCMech)

#else

#ifndef MPM_FIX_NEUMANN_BC_MECH_H
#define MPM_FIX_NEUMANN_BC_MECH_H

#include <fix.h>
#include <var.h>
#include <vector>

class FixNeumannBCMech : public Fix {
 public:
  FixNeumannBCMech(class MPM *, vector<string>);
  ~FixNeumannBCMech();
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

  void write_restart(ofstream *);
  void read_restart(ifstream *);

private:
  const map<int, string> usage = {
      {1, "Usage: fix(fix-ID, neumann_bc_mech, group, tx)\n"},
      {2, "Usage: fix(fix-ID, neumann_bc_mech, group, tx, ty)\n"},
      {3, "Usage: fix(fix-ID, neumann_bc_mech, group, tx, ty, tz)\n"}};
  const map<int, int> Nargs = {{1, 4}, {2, 5}, {3, 6}};

  class Var t[3];                     //< Surface traction/compression pressure
};

#endif
#endif

