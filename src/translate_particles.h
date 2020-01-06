/* -*- c++ -*- ----------------------------------------------------------*/

#ifdef COMMAND_CLASS

CommandStyle(translate_particles,TranslateParticles)

#else

#ifndef MPM_TRANSLATE_PARTICLES_H
#define MPM_TRANSLATE_PARTICLES_H

#include "var.h"
#include "pointers.h"

class TranslateParticles : protected Pointers {
 public:
  TranslateParticles(class MPM *);
  class Var command(vector<string>);

 private:
  void translate_region(vector<string>, int);
  class Var delx, dely, delz;    // Set displacements in x, y, and z directions.
  bool xset, yset, zset;
};

#endif
#endif
