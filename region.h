/* -*- c++ -*- ----------------------------------------------------------*/

#ifndef MPM_REGION_H
#define MPM_REGION_H

#include "pointers.h"

class Region : protected Pointers {
 public:
  Region(class MPM *, string *);
  virtual ~Region();

  //protected:
};

#endif
