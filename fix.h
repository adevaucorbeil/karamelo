/* -*- c++ -*- ----------------------------------------------------------*/

#ifndef MPM_FIX_H
#define MPM_FIX_H

#include "pointers.h"
#include <vector>

class Fix : protected Pointers {
 public:
  string id;

  Fix(class MPM *, vector<string>);
  virtual ~Fix() {};
  virtual void init() = 0;
  virtual void setup() = 0;
};

#endif
