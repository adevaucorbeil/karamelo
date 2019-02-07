/* -*- c++ -*- ----------------------------------------------------------*/

#ifndef MPM_COMPUTE_H
#define MPM_COMPUTE_H

#include "pointers.h"
#include <vector>

class Compute : protected Pointers {
 public:
  string id;
  int igroup;

  Compute(class MPM *, vector<string>);
  virtual ~Compute() {};
  virtual void init() = 0;
  virtual void setup() = 0;
};

#endif
