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

#ifndef MPM_DAMAGE_H
#define MPM_DAMAGE_H

#include <pointers.h>
#include <vector>
#include <matrix.h>
#include <grid.h>

class Solid;

class Damage : protected Pointers {
 public:
  string id;                        ///< damage identification string
  string style;                     ///< damage style

  Damage(MPM *, vector<string>);
  virtual ~Damage();
  virtual void init();
  void options(vector<string> *, vector<string>::iterator);

  virtual void write_restart(ofstream*) = 0;
  virtual void read_restart(ifstream*) = 0;

  // implemented by each Damage
  virtual void compute_damage(Solid &solid,
                              Kokkos::View<double*> &pH,
                              Kokkos::View<Matrix3d*> &sigma_dev,
                              Kokkos::View<double*> &plastic_strain_increment) const = 0;
};

#endif
