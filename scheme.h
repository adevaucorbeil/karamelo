/* -*- c++ -*- ----------------------------------------------------------*/

#ifndef MPM_SCHEME_H
#define MPM_SCHEME_H

#include "pointers.h"

class Scheme : protected Pointers {
 public:
  Scheme(class MPM *);
  virtual ~Scheme();
  virtual void setup() = 0;
  virtual void run(class Var) = 0;
};

#endif
