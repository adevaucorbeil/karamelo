/* -*- c++ -*- ----------------------------------------------------------*/

#ifndef MPM_SOLID_H
#define MPM_SOLID_H

#include "pointers.h"

class Solid : protected Pointers {
 public:
  int np;

  // functions
  Solid(class MPM *);
  ~Solid();

};

#endif
