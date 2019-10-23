/* -*- c++ -*- ----------------------------------------------------------*/
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
  protected:
  vector<string> known_var = {"x", "y", "z",
			      "x0", "y0", "z0",
			      "vx", "vy", "vz",
			      "s11", "s22", "s33",
			      "s12", "s13", "s23",
			      "seq", "volume", "mass",
			      "damage", "damage_init",
			      "bx", "by", "bz",
			      "ep", "epdot", "T", "ienergy"};
};

#endif
#endif
