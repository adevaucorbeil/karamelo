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

FixStyle(initial_velocity_nodes,FixInitialVelocityNodes)

#else

#ifndef MPM_FIX_INITIAL_VELOCITY_NODES_H
#define MPM_FIX_INITIAL_VELOCITY_NODES_H

#include <fix.h>
#include <matrix.h>

class Expression;

class FixInitialVelocityNodes : public Fix {
 public:
  FixInitialVelocityNodes(MPM *, vector<string>);

  void prepare();

  void post_update_grid_state(Grid &grid);
  void post_velocities_to_grid(Grid &grid);

  void write_restart(ofstream *) {};
  void read_restart(ifstream *) {};

private:
  const map<int, string> usage = {
      {1, "Usage: fix(fix-ID, initial_velocity_nodes, group, vx)\n"},
      {2, "Usage: fix(fix-ID, initial_velocity_nodes, group, vx, vy)\n"},
      {3, "Usage: fix(fix-ID, initial_velocity_nodes, group, vx, vy, vz)\n"}};
  const map<int, int> Nargs = {{1, 4}, {2, 5}, {3, 6}};

  Expression *v[3];
};

#endif
#endif

