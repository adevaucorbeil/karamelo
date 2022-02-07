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

  void setup(vector<string>);

  void compute_grid_weight_functions_and_gradients();
  double (*basis_function)(double, int);
  double (*derivative_basis_function)(double, int, double);

  vector<Grid *> grids() override;
  bool should_compute_mass_nodes() override;
  void compute_internal_force_nodes(Solid &solid, int in, int ip, double wf, const Vector3d &wfd) override;
  void check_particle_in_domain(const Vector3d &x, int ip) override;

  void compute_rate_deformation_gradient(bool);
  void update_deformation_gradient();
  void update_stress(bool);
  void adjust_dt();
  void reset();
  void exchange_particles();

private:
  int update_Di;
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
