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

  float G();

  void write_restart(ofstream *);
  void read_restart(ifstream *);

  void update_deviatoric_stress(Solid &solid,
                                Kokkos::View<float*> &plastic_strain_increment,
                                Kokkos::View<Matrix3d*> &sigma_dev) const override;

protected:
  float G_, yieldStress;
};

#endif
#endif
