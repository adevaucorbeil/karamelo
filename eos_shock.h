/* -*- c++ -*- ----------------------------------------------------------*/


#ifdef EOS_CLASS

EOSStyle(shock,EOSShock)

#else

#ifndef MPM_EOS_SHOCK_H
#define MPM_EOS_SHOCK_H

#include "eos.h"
#include <Eigen/Eigen>

class EOSShock : public EOS {

public:
  EOSShock(class MPM *, vector<string>);
  ~EOSShock();

  double rho0();
  double K();
  double G();
  double compute_pressure(const double, const double, const double);

protected:
  double rho0_, K_, e0, c0, S, Gamma;
};

#endif
#endif
