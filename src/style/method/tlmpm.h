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

MethodStyle(tlmpm,TLMPM)

#else

#ifndef LMP_TLMPM_H
#define LMP_TLMPM_H

#include <method.h>
#include <vector>
#include <matrix.h>

class TLMPM : public Method {
public:
  TLMPM(class MPM *);

  void setup(vector<string>);
  void compute_grid_weight_functions_and_gradients(Solid &solid, int ip) override;
  vector<Grid *> grids() override;
  void reset_mass_nodes(Grid &grid, int in) override;
  void compute_mass_nodes(Solid &solid, int ip) override;
  void compute_internal_force_nodes(Solid &solid, int ip) override;
  void update_grid_positions(Grid &grid, int in) override;
  vector<Matrix3d> &get_gradients(Solid &solid) override;
  virtual void update_deformation_gradient_matrix(Solid &solid, int ip) override;
  void update_velocity_gradient_matrix(Solid &solid, int ip) override;
};

// double linear_basis_function(double, int);
// double derivative_linear_basis_function(double, int, double);
// double cubic_spline_basis_function(double, int);
// double derivative_cubic_spline_basis_function(double, int, double);
// double bernstein_quadratic_basis_function(double, int);
// double derivative_bernstein_quadratic_basis_function(double, int, double);


#endif
#endif
