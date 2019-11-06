/* -*- c++ -*- ----------------------------------------------------------*/


#ifdef STRENGTH_CLASS

StrengthStyle(plastic,StrengthPlastic)

#else

#ifndef MPM_STRENGTH_PLASTIC_H
#define MPM_STRENGTH_PLASTIC_H

#include "strength.h"
#include <Eigen/Eigen>

class StrengthPlastic : public Strength {

public:
  StrengthPlastic(class MPM *, vector<string>);
  ~StrengthPlastic() {};

  double G();

  Eigen::Matrix3d  update_deviatoric_stress
  ( const Eigen::Matrix3d& sigma,
    const Eigen::Matrix3d& D,
    double &               plastic_strain_increment,
    const double           eff_plastic_strain,
    const double           epsdot,
    const double           damage,
    const double           temperature);

protected:
  double G_, yieldStress;
};

#endif
#endif
