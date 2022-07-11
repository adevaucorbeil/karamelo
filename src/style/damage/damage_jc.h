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


#ifdef DAMAGE_CLASS

DamageStyle(damage_johnson_cook,DamageJohnsonCook)

#else

#ifndef MPM_DAMAGE_JOHNSON_COOK_H
#define MPM_DAMAGE_JOHNSON_COOK_H

#include <damage.h>
#include <matrix.h>

class DamageJohnsonCook : public Damage {

public:
  DamageJohnsonCook(class MPM *, vector<string>);
  ~DamageJohnsonCook() {};

  void write_restart(ofstream *);
  void read_restart(ifstream *);

  void compute_damage(Solid &solid,
                      Kokkos::View<float*> &pH,
                      Kokkos::View<Matrix3d*> &sigma_dev,
                      Kokkos::View<float*> &plastic_strain_increment) const override;

protected:
  float d1, d2, d3, d4, d5, epsdot0, Tr, Tm, Tmr;
  string usage = "Usage: damage(damage-ID, damage_johnson_cook, d1, d2, d3, d4, d5, epsdot0, Tr, Tm)\n";
  int Nargs = 10;
};

#endif
#endif
