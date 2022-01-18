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

FixStyle(convection_particles,FixConvectionParticles)

#else

#ifndef MPM_FIX_CONVECTION_PARTICLES_H
#define MPM_FIX_CONVECTION_PARTICLES_H

#include "fix.h"
#include "var.h"
#include <Eigen/Eigen>
#include <vector>

class FixConvectionParticles : public Fix {
 public:
  FixConvectionParticles(class MPM *, vector<string>);
  ~FixConvectionParticles();
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
  string usage = "Usage: fix(fix-ID, convection_particles, group, h, Tinf)\n";
  const int Nargs = 5;

  class Var Tinf;                  //< Ambiant temperature.
  double h;                        //< Heat transfer coefficient.


};

#endif
#endif

