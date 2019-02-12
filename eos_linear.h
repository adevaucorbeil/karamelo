/* -*- c++ -*- ----------------------------------------------------------*/


#ifdef EOS_CLASS

EOSStyle(linear,EOSLinear)

#else

#ifndef MPM_EOS_BLOCK_H
#define MPM_EOS_BLOCK_H

#include "eos.h"
#include <Eigen/Eigen>

class EOSLinear : public EOS {

public:
  EOSLinear(class MPM *, vector<string>);
  ~EOSLinear();

  double rho0();
  double K();
  double G();
  double compute_pressure(double);
  Eigen::Matrix3d update_deviatoric_stress(Eigen::Matrix3d, Eigen::Matrix3d);
  void update_stress(Eigen::Matrix3d&, Eigen::Matrix3d, double);

protected:
  double rho0_, K_, G_;
};

#endif
#endif
