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

FixStyle(velocity_nodes,FixVelocityNodes)

#else

#ifndef MPM_FIX_VELOCITY_NODES_H
#define MPM_FIX_VELOCITY_NODES_H

#include "fix.h"
#include "var.h"
#include <vector>

class FixVelocityNodes : public Fix {
 public:
  FixVelocityNodes(class MPM *, vector<string>);
  ~FixVelocityNodes();
  void setmask();
  void init();
  void setup();
  
  void initial_integrate() {};
  void post_particles_to_grid() {};
  void post_update_grid_state();
  void post_grid_to_point() {};
  void post_advance_particles() {};
  void post_velocities_to_grid();
  void final_integrate() {};

  void write_restart(ofstream *);
  void read_restart(ifstream *);

private:
  const map<int, string> usage = {
      {1, "Usage: fix(fix-ID, velocity_nodes, group, vx)\n"},
      {2, "Usage: fix(fix-ID, velocity_nodes, group, vx, vy)\n"},
      {3, "Usage: fix(fix-ID, velocity_nodes, group, vx, vy, vz)\n"}};
  const map<int, int> Nargs = {{1, 4}, {2, 5}, {3, 6}};

  class Var xvalue, yvalue, zvalue;                  //< Velocities in x, y, and z directions.
  class Var xprevvalue, yprevvalue, zprevvalue;      //< Velocities in x, y, and z directions from previous time step.
  bool xset, yset, zset;                             //< Does the fix set the x, y, and z velocities of the group?
};

#endif
#endif

