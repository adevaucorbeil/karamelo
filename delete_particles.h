/* -*- c++ -*- ----------------------------------------------------------*/

#ifdef COMMAND_CLASS

CommandStyle(delete_particles,DeleteParticles)

#else

#ifndef MPM_DELETE_PARTICLES_H
#define MPM_DELETE_PARTICLES_H

#include "var.h"
#include "pointers.h"

class DeleteParticles : protected Pointers {
 public:
  DeleteParticles(class MPM *);
  class Var command(vector<string>);

 private:
  int **dlist;
  void delete_region(vector<string>, int);
};

#endif
#endif
