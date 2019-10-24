/* -*- c++ -*- ----------------------------------------------------------*/


#ifdef STRENGTH_CLASS

StrengthStyle(fluid,StrengthFluid)

#else

#ifndef MPM_STRENGTH_FLUID_H
#define MPM_STRENGTH_FLUID_H

#include "strength.h"
#include <Eigen/Eigen>

class StrengthFluid : public Strength {

public:
  StrengthFluid(class MPM *, vector<string>);
  ~StrengthFluid() {};

  double G();

  Eigen::Matrix3d update_deviatoric_stress
    (const Eigen::Matrix3d& sigma,
     const Eigen::Matrix3d& D,
     double &plastic_strain_increment,
     const double eff_plastic_strain,
     const double epsdot,
     const double damage);

protected:
  double G_;
};

#endif
#endif
