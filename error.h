/* -*- c++ -*- ---------------------------------------------------------- */

#ifndef MPM_ERROR_H
#define MPM_ERROR_H

#include "pointers.h"

class Error : protected Pointers {
 public:
  Error(class MPM *);

  void all(const char *, int, const string);
  void one(const char *, int, const string);
};

#endif
