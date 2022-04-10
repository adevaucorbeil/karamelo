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

#ifndef MPM_STRENGTH_H
#define MPM_STRENGTH_H

#include <pointers.h>
#include <vector>
#include <matrix.h>
#include <grid.h>

class Solid;

class Strength : protected Pointers {
 public:
  string id;                        ///< Strength identification string
  string style;                     ///< Strength style

  Strength(MPM *, vector<string>);
  virtual ~Strength();
  virtual void init();
  void options(vector<string> *, vector<string>::iterator);

  virtual void write_restart(ofstream*) = 0;
  virtual void read_restart(ifstream*) = 0;

  // implemented by each Strength
  //virtual compute_pressure()
  virtual double G() = 0;

  virtual void update_deviatoric_stress(Solid &solid,
                                        Kokkos::View<double*, MemorySpace> &plastic_strain_increment,
                                        Kokkos::View<Matrix3d*, MemorySpace> &sigma_dev) const = 0;
  //protected:
};

#endif
