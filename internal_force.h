/* -*- c++ -*- ----------------------------------------------------------*/

#ifdef COMMAND_CLASS

CommandStyle(internal_force,IntForce)

#else

#ifndef LMP_INTERNAL_FORCE_H
#define LMP_INTERNAL_FORCE_H

#include "pointers.h"

class IntForce : protected Pointers {
 public:
  IntForce(class MPM *);
  class Var command(vector<string>);

 private:
  int igroup, dir;
};

#endif
#endif
