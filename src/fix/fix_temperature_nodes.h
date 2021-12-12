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
#include <vector>

class FixTemperatureNodes : public Fix {
 public:
  FixTemperatureNodes(class MPM *, vector<string>);
  ~FixTemperatureNodes();
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
  string usage = "Usage: fix(fix-ID, temperature_nodes, group-ID, T)\n";
  int Nargs = 4;

  class Var Tvalue;                      //< Temperature variable.
  class Var Tprevvalue;                  //< Temperature variable from previous time step.
};

#endif
#endif

