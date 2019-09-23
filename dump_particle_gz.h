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
  //protected:
};

#endif
#endif
