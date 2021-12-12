/* -*- c++ -*- ----------------------------------------------------------*/

#ifdef COMMAND_CLASS

CommandStyle(translate_particles,TranslateParticles)

#else

#ifndef MPM_TRANSLATE_PARTICLES_H
#define MPM_TRANSLATE_PARTICLES_H

#include <pointers.h>
#include <var.h>

class TranslateParticles : protected Pointers {
 public:
  TranslateParticles(class MPM *);
  class Var command(vector<string>);

 private:
  string usage = "Usage: translate_particles(solid-ID, region, region-ID, delx, dely, delz)\n";
  int Nargs = 6;
  void translate_region(vector<string>, int);
  class Var delx, dely, delz;    // Set displacements in x, y, and z directions.
  bool xset, yset, zset;
};

#endif
#endif
