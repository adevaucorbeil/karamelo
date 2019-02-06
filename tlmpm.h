/* -*- c++ -*- ----------------------------------------------------------*/

#ifdef METHOD_CLASS

MethodStyle(tlmpm,TLMPM)

#else

#ifndef LMP_TLMPM_H
#define LMP_TLMPM_H

#include "method.h"
#include <vector>

class TLMPM : public Method {
 public:
  TLMPM(class MPM *, vector<string>);
  ~TLMPM() {}

 protected:
};

#endif
#endif
