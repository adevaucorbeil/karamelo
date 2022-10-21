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

#ifdef DUMP_CLASS

DumpStyle(particle/gz,DumpParticleGz)

#else

#ifndef MPM_DUMP_PARTICLE_GZ_H
#define MPM_DUMP_PARTICLE_GZ_H

#include "dump.h"

class DumpParticleGz : public Dump {
 public:
  DumpParticleGz(MPM *, vector<string>);
  ~DumpParticleGz();

  void write();
  protected:
  vector<string> known_var = {"x", "y", "z",
			      "x0", "y0", "z0",
			      "vx", "vy", "vz",
			      "s11", "s22", "s33",
			      "s12", "s13", "s23",
			      "e11", "e22", "e33",
			      "e12", "e13", "e23",
			      "seq", "volume", "mass",
			      "damage", "damage_init",
			      "bx", "by", "bz",
			      "ep", "epdot", "T",
			      "ienergy", "gamma", "surf"};
};

#endif
#endif
