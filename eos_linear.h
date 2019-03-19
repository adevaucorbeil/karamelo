/* -*- c++ -*- ----------------------------------------------------------*/


#ifdef EOS_CLASS

EOSStyle(linear,EOSLinear)

#else

#ifndef MPM_EOS_LINEAR_H
#define MPM_EOS_LINEAR_H

#include "eos.h"
#include <Eigen/Eigen>

class EOSLinear : public EOS {

public:
  EOSLinear(class MPM *, vector<string>);
  ~EOSLinear();

  double rho0();
  double K();
  double G();
  double compute_pressure(const double J, const double rho, const double e);

protected:
  double rho0_, K_;
};

#endif
#endif
