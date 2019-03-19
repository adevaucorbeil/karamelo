/* -*- c++ -*- ----------------------------------------------------------*/


#ifdef STRENGTH_CLASS

StrengthStyle(linear,StrengthLinear)

#else

#ifndef MPM_STRENGTH_LINEAR_H
#define MPM_STRENGTH_LINEAR_H

#include "strength.h"
#include <Eigen/Eigen>

class StrengthLinear : public Strength {

public:
  StrengthLinear(class MPM *, vector<string>);
  ~StrengthLinear() {};

  double G();
  Eigen::Matrix3d update_deviatoric_stress(const Eigen::Matrix3d, const Eigen::Matrix3d);

protected:
  double G_;
};

#endif
#endif
