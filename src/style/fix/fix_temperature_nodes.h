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

FixStyle(temperature_nodes,FixTemperatureNodes)

#else

#ifndef MPM_FIX_TEMPERATURE_NODES_H
#define MPM_FIX_TEMPERATURE_NODES_H

#include <fix.h>
#include <var.h>
#include <matrix.h>

class FixTemperatureNodes : public Fix {
 public:
  FixTemperatureNodes(MPM *, vector<string>);

  void prepare();
  void reduce();
  
  void post_update_grid_state();
  void post_velocities_to_grid();

  void write_restart(ofstream *);
  void read_restart(ifstream *);

private:  
  string usage = "Usage: fix(fix-ID, temperature_nodes, group-ID, T)\n";
  int Nargs = 4;
  string Targ;

  Var Tvalue;                      //< Temperature variable.
  Var Tprevvalue;                  //< Temperature variable from previous time step.
};

#endif
#endif

