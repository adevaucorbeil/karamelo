/* -*- c++ -*- ----------------------------------------------------------
 *
 *                    ***       Karamelo       ***
 *               Parallel Material Point Method Simulator
 * 
 * Copyright (2019) Alban de Vaucorbeil, alban.devaucorbeil@monash.edu
 * Materials Science and Engineering, Monash University
 * Clayton VIC 3800, Australia

 * This software is distributed under the GNU General Public License.
 *
 * ----------------------------------------------------------------------- */


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
  
  Eigen::Matrix3d  update_deviatoric_stress
  ( const Eigen::Matrix3d& sigma,
    const Eigen::Matrix3d& D,
    double &               plastic_strain_increment,
    const double           eff_plastic_strain,
    const double           epsdot,
    const double           damage,
    const double           temperature);

protected:
  double G_, A, B, n, m, epsdot0, C, Tr, Tm, Tmr;
  string usage = "Usage: strength(strength-ID, johnson_cook, G, A, B, n, epsdot0, C, m, Tr, Tm)\n";
  int Nargs = 11;
};

#endif
#endif
