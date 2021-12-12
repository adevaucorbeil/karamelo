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

StrengthStyle(plastic,StrengthPlastic)

#else

#ifndef MPM_STRENGTH_PLASTIC_H
#define MPM_STRENGTH_PLASTIC_H

#include <strength.h>
#include <matrix.h>

class StrengthPlastic : public Strength {

public:
  StrengthPlastic(class MPM *, vector<string>);
  ~StrengthPlastic() {};

  double G();

  void write_restart(ofstream *);
  void read_restart(ifstream *);

  Matrix3d  update_deviatoric_stress
  ( const Matrix3d& sigma,
    const Matrix3d& D,
    double &               plastic_strain_increment,
    const double           eff_plastic_strain,
    const double           epsdot,
    const double           damage,
    const double           temperature = 0);

protected:
  double G_, yieldStress;
};

#endif
#endif