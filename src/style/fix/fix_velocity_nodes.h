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

#include <fix.h>
#include <var.h>
#include <matrix.h>

class FixVelocityNodes : public Fix {
 public:
  FixVelocityNodes(MPM *, vector<string>);

  void prepare();
  void reduce();
  
  void post_update_grid_state(Grid &grid, int in);
  void post_velocities_to_grid(Grid &grid, int in);

  void write_restart(ofstream *);
  void read_restart(ifstream *);

private:
  const map<int, string> usage = {
      {1, "Usage: fix(fix-ID, velocity_nodes, group, vx)\n"},
      {2, "Usage: fix(fix-ID, velocity_nodes, group, vx, vy)\n"},
      {3, "Usage: fix(fix-ID, velocity_nodes, group, vx, vy, vz)\n"}};
  const map<int, int> Nargs = {{1, 4}, {2, 5}, {3, 6}};

  Var xvalue, yvalue, zvalue;                  //< Velocities in x, y, and z directions.
  Var xprevvalue, yprevvalue, zprevvalue;      //< Velocities in x, y, and z directions from previous time step.
  bool xset, yset, zset;                             //< Does the fix set the x, y, and z velocities of the group?
  Vector3d ftot;
};

#endif
#endif

