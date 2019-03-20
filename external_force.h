/* -*- c++ -*- ----------------------------------------------------------*/

#ifdef COMMAND_CLASS

CommandStyle(external_force,ExtForce)

#else

#ifndef LMP_EXTERNAL_FORCE_H
#define LMP_EXTERNAL_FORCE_H

#include "pointers.h"

class ExtForce : protected Pointers {
 public:
  ExtForce(class MPM *);
  class Var command(vector<string>);

 private:
  int igroup, dir;
};

#endif
#endif
