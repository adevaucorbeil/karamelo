/* -*- c++ -*- ----------------------------------------------------------*/

#ifdef COMMAND_CLASS

CommandStyle(xcm,CentreOfMass)

#else

#ifndef LMP_CENTRE_OF_MASS_H
#define LMP_CENTRE_OF_MASS_H

#include "pointers.h"

class CentreOfMass : protected Pointers {
 public:
  CentreOfMass(class MPM *);
  class Var command(vector<string>);

 private:
  int igroup, dir;
};

#endif
#endif
