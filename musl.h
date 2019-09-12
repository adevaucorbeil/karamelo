/* -*- c++ -*- ----------------------------------------------------------*/

#ifdef SCHEME_CLASS

SchemeStyle(musl,MUSL)

#else

#ifndef LMP_MUSL_H
#define LMP_MUSL_H

#include "scheme.h"
#include "var.h"
#include <vector>

class MUSL : public Scheme {
 public:
  MUSL(class MPM *, vector<string>);
  ~MUSL() {}
  void setup();
  void run(class Var);
};

#endif
#endif
