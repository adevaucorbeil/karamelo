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

#include "method.h"
#include <vector>
#include <Eigen/Eigen>


class ULMPM : public Method {
 public:
  bool apic;
  
  ULMPM(class MPM *);
  ~ULMPM();

  void setup(vector<string>);

  void compute_grid_weight_functions_and_gradients();
  double (*basis_function)(double, int);
  double (*derivative_basis_function)(double, int, double);
  void particles_to_grid();
  void update_grid_state();
  void grid_to_points();
  void advance_particles();
  void velocities_to_grid();
  void update_grid_positions() {};
  void compute_rate_deformation_gradient();
  void update_deformation_gradient();
  void update_stress();
  void adjust_dt();
  void reset();
  void exchange_particles();

private:
  int update_wf;
  int rigid_solids;
};

// double linear_basis_function(double, int);
// double derivative_linear_basis_function(double, int, double);
// double cubic_spline_basis_function(double, int);
// double derivative_cubic_spline_basis_function(double, int, double);
// double bernstein_quadratic_basis_function(double, int);
// double derivative_bernstein_quadratic_basis_function(double, int, double);


#endif
#endif
