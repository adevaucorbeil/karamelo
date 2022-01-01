/* -*- c++ -*- ----------------------------------------------------------*/

#ifdef METHOD_CLASS

MethodStyle(tlcpdi,TLCPDI)

#else

#ifndef LMP_TLCPDI_H
#define LMP_TLCPDI_H

#include <method.h>
#include <vector>
#include <matrix.h>


class TLCPDI : public Method {
 public:

  TLCPDI(class MPM *);
  ~TLCPDI();

  void setup(vector<string>);

  void compute_grid_weight_functions_and_gradients();
  double (*basis_function)(double, int);
  double (*derivative_basis_function)(double, int, double);
  void particles_to_grid();
  void particles_to_grid_USF_1();
  void particles_to_grid_USF_2();
  void update_grid_state();
  void grid_to_points();
  void advance_particles();
  void velocities_to_grid();
  void update_grid_positions();
  void compute_rate_deformation_gradient(bool);
  void update_deformation_gradient();
  void update_stress(bool);
  void adjust_dt();
  void reset();
  void exchange_particles() {};

  bool update_wf;

  private:
  string known_styles[2] = {"R4", "Q4"};
};

// double linear_basis_function(double, int);
// double derivative_linear_basis_function(double, int, double);
// double cubic_spline_basis_function(double, int);
// double derivative_cubic_spline_basis_function(double, int, double);
// double bernstein_quadratic_basis_function(double, int);
// double derivative_bernstein_quadratic_basis_function(double, int, double);


#endif
#endif
