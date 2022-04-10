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

StrengthStyle(linear,StrengthLinear)

#else

#ifndef MPM_STRENGTH_LINEAR_H
#define MPM_STRENGTH_LINEAR_H

#include <strength.h>
#include <matrix.h>

class StrengthLinear : public Strength {

public:
  StrengthLinear(class MPM *, vector<string>);
  ~StrengthLinear() {};

  double G();

  void write_restart(ofstream *);
  void read_restart(ifstream *);

  void update_deviatoric_stress(Solid &solid,
                                Kokkos::View<double*, MemorySpace> &plastic_strain_increment,
                                Kokkos::View<Matrix3d*, MemorySpace> &sigma_dev) const override;

protected:
  double G_;
};

#endif
#endif
