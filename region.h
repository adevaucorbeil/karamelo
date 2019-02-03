/* -*- c++ -*- ----------------------------------------------------------*/

#ifndef MPM_REGION_H
#define MPM_REGION_H

#include "pointers.h"
#include <vector>

class Region : protected Pointers {
 public:
  string id;

  Region(class MPM *, vector<string>);
  virtual ~Region();
  virtual void init();

  //protected:
};

#endif
