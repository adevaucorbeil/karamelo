/* -*- c++ -*- ----------------------------------------------------------*/

#ifndef MPM_METHOD_H
#define MPM_METHOD_H

#include "pointers.h"

class Method : protected Pointers {
 public:
  Method(class MPM *);
  virtual ~Method();
};

#endif
