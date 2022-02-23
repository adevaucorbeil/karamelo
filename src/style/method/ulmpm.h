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

#ifdef METHOD_CLASS

MethodStyle(ulmpm,ULMPM)

#else

#ifndef LMP_ULMPM_H
#define LMP_ULMPM_H

#include <method.h>
#include <vector>
#include <matrix.h>


class ULMPM : public Method {
 public:
  bool apic;
  
  ULMPM(class MPM *);

  void setup(vector<string> args) override;
  vector<Grid *> grids() override;
  void compute_internal_force_nodes(Solid &solid, int ip) override;
  void check_particle_in_domain(const Vector3d &x, int ip) override;
  Kokkos::View<Matrix3d*> &get_gradients(Solid &solid) override;
  virtual void update_deformation_gradient_matrix(Solid &solid, int ip) override;
  void update_velocity_gradient_matrix(Solid &solid, int ip) override;
  void exchange_particles() override;
};

// double linear_basis_function(double, int);
// double derivative_linear_basis_function(double, int, double);
// double cubic_spline_basis_function(double, int);
// double derivative_cubic_spline_basis_function(double, int, double);
// double bernstein_quadratic_basis_function(double, int);
// double derivative_bernstein_quadratic_basis_function(double, int, double);


#endif
#endif
