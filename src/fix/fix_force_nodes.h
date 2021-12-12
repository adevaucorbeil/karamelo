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

FixStyle(force_nodes,FixForceNodes)

#else

#ifndef MPM_FIX_FORCE_NODES_H
#define MPM_FIX_FORCE_NODES_H

#include <fix.h>
#include <var.h>
#include <vector>

class FixForceNodes : public Fix {
 public:
  FixForceNodes(class MPM *, vector<string>);
  ~FixForceNodes();
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

  void write_restart(ofstream *);
  void read_restart(ifstream *);

private:
  class Var xvalue, yvalue, zvalue;    // Set force in x, y, and z directions.
  bool xset, yset, zset;               // Does the fix set the x, y, and z forces of the group?
};

#endif
#endif

