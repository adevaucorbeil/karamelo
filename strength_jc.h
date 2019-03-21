/* -*- c++ -*- ----------------------------------------------------------*/


#ifdef STRENGTH_CLASS

StrengthStyle(johnson_cook,StrengthJohnsonCook)

#else

#ifndef MPM_STRENGTH_JOHNSON_COOK_H
#define MPM_STRENGTH_JOHNSON_COOK_H

#include "strength.h"
#include <Eigen/Eigen>

class StrengthJohnsonCook : public Strength {

public:
  StrengthJohnsonCook(class MPM *, vector<string>);
  ~StrengthJohnsonCook() {};

  double G();
  Eigen::Matrix3d update_deviatoric_stress(const Eigen::Matrix3d sigma, const Eigen::Matrix3d D, double &plastic_strain_increment, const double eff_plastic_strain, const double epsdot);

protected:
  double G_, A, B, n, epsdot0, C;
};

#endif
#endif
