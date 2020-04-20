/* -*- c++ -*- ----------------------------------------------------------*/


#ifdef EOS_CLASS

EOSStyle(fluid,EOSFluid)

#else

#ifndef MPM_EOS_FLUID_H
#define MPM_EOS_FLUID_H

#include "eos.h"
#include <Eigen/Eigen>

class EOSFluid : public EOS {

public:
  EOSFluid(class MPM *, vector<string>);
  ~EOSFluid();

  double rho0();
  double K();
  double G();
  void compute_pressure(double &, double &, const double, const double, const double, const double, const Eigen::Matrix3d, const double);

protected:
  double rho0_, K_, Gamma;
};

#endif
#endif
