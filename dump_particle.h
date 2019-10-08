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

DumpStyle(particle,DumpParticle)

#else

#ifndef MPM_DUMP_PARTICLE_H
#define MPM_DUMP_PARTICLE_H

#include "dump.h"

class DumpParticle : public Dump {
 public:
  DumpParticle(MPM *, vector<string>);
  ~DumpParticle();

  void write();
  //protected:
};

#endif
#endif
