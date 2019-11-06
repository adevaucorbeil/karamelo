/* -*- c++ -*- ----------------------------------------------------------*/

#ifndef MPM_STRENGTH_H
#define MPM_STRENGTH_H

#include "pointers.h"
#include <vector>
#include <Eigen/Eigen>

class Strength : protected Pointers {
 public:
  string id;

  Strength(MPM *, vector<string>);
  virtual ~Strength();
  virtual void init();
  void options(vector<string> *, vector<string>::iterator);

  // implemented by each Strength
  //virtual compute_pressure()
  virtual double G() = 0;

  virtual Eigen::Matrix3d  update_deviatoric_stress
  ( const Eigen::Matrix3d& sigma,
    const Eigen::Matrix3d& D,
    double &               plastic_strain_increment,
    const double           eff_plastic_strain,
    const double           epsdot,
    const double           damage,
    const double           temperature) = 0;
  //protected:
};

#endif
