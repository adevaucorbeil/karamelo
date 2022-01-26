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

#include <fix.h>
#include <var.h>

class FixConvectionParticles : public Fix {
 public:
  FixConvectionParticles(MPM *, vector<string>);

  void prepare();
  void reduce();
  
  void initial_integrate();

  void write_restart(ofstream *);
  void read_restart(ifstream *);

private:
  string usage = "Usage: fix(fix-ID, convection_particles, group, h, Tinf)\n";
  const int Nargs = 5;

  Var Tinf;                  //< Ambiant temperature.
  double h;                        //< Heat transfer coefficient.
  double qtot;
};

#endif
#endif

