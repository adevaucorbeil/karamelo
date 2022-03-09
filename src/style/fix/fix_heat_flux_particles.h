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

FixStyle(heat_flux_particles,FixHeatFluxParticles)

#else

#ifndef MPM_FIX_HEAT_FLUX_PARTICLES_H
#define MPM_FIX_HEAT_FLUX_PARTICLES_H

#include <fix.h>
#include <var.h>

class FixHeatFluxParticles : public Fix {
 public:
  FixHeatFluxParticles(MPM *, vector<string>);

  void prepare();
  void reduce();
  
  void initial_integrate(Solid &solid);

  void write_restart(ofstream *);
  void read_restart(ifstream *);

private:
  string usage = "Usage: fix(fix-ID, heat_flux_particles, group, q)\n";
  const int Nargs = 4;

  Var q;                  //< Flux
  double qtot;
};

#endif
#endif

