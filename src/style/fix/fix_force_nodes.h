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
#include <matrix.h>

class FixForceNodes : public Fix {
 public:
  FixForceNodes(MPM *, vector<string>);

  void prepare();
  void reduce();
  
  void post_particles_to_grid();

  void write_restart(ofstream *);
  void read_restart(ifstream *);

private:
  Var xvalue, yvalue, zvalue;    // Set force in x, y, and z directions.
  bool xset, yset, zset;               // Does the fix set the x, y, and z forces of the group?
  Vector3d ftot;
};

#endif
#endif

