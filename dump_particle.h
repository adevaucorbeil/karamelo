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

  //protected:
};

#endif
#endif
