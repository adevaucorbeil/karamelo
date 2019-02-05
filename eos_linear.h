/* -*- c++ -*- ----------------------------------------------------------*/


#ifdef EOS_CLASS

EOSStyle(linear,EOSLinear)

#else

#ifndef MPM_EOS_BLOCK_H
#define MPM_EOS_BLOCK_H

#include "eos.h"
#include <stdarg.h>

class EOSLinear : public EOS {

 public:
  EOSLinear(class MPM *, vector<string>);
  ~EOSLinear();

  double rho0();
  double K();
  double compute_pressure(double);

 protected:
  double rho0_, K_;
};

#endif
#endif
