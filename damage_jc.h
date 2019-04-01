/* -*- c++ -*- ----------------------------------------------------------*/


#ifdef DAMAGE_CLASS

DamageStyle(damage_johnson_cook,DamageJohnsonCook)

#else

#ifndef MPM_DAMAGE_JOHNSON_COOK_H
#define MPM_DAMAGE_JOHNSON_COOK_H

#include "damage.h"
#include <Eigen/Eigen>

class DamageJohnsonCook : public Damage {

public:
  DamageJohnsonCook(class MPM *, vector<string>);
  ~DamageJohnsonCook() {};

  void compute_damage(double &damage_init, double &damage, const double pH, const Eigen::Matrix3d Sdev, const double epsdot, const double plastic_strain_increment);

protected:
  double d1, d2, d3, d4, epsdot0;
};

#endif
#endif
