/* -*- c++ -*- ----------------------------------------------------------*/

#ifndef MPM_SCHEME_H
#define MPM_SCHEME_H

#include "pointers.h"

class Scheme : protected Pointers {
 public:
  Scheme(class MPM *);
  virtual ~Scheme();
  virtual void run(int) = 0;
};

#endif
