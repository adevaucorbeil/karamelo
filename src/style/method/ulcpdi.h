///* -*- c++ -*- ----------------------------------------------------------
// *
// *                    ***       Karamelo       ***
// *               Parallel Material Point Method Simulator
// * 
// * Copyright (2019) Alban de Vaucorbeil, alban.devaucorbeil@monash.edu
// * Materials Science and Engineering, Monash University
// * Clayton VIC 3800, Australia
//
// * This software is distributed under the GNU General Public License.
// *
// * ----------------------------------------------------------------------- */
//
//#ifdef METHOD_CLASS
//
//MethodStyle(ulcpdi,ULCPDI)
//
//#else
//
//#ifndef LMP_ULCPDI_H
//#define LMP_ULCPDI_H
//
//#include <method.h>
//#include <vector>
//#include <matrix.h>
//#include <map>
//
//
//class ULCPDI : public Method {
// public:
//  float FLIP;
//  string shape_function;
//
//  ULCPDI(class MPM *);
//  ~ULCPDI();
//
//  void setup(vector<string>);
//
//  void compute_grid_weight_functions_and_gradients();
//  float (*basis_function)(float, int);
//  float (*derivative_basis_function)(float, int, float);
//  void particles_to_grid();
//  void particles_to_grid_USF_1();
//  void particles_to_grid_USF_2();
//  void update_grid_state();
//  void grid_to_points();
//  void advance_particles();
//  void velocities_to_grid();
//  void update_grid_positions() {};
//  void compute_rate_deformation_gradient(bool);
//  void update_deformation_gradient();
//  void update_stress(bool);
//  void adjust_dt();
//  void reset();
//  void exchange_particles();
//
//  int update_wf;
//private:
//  string known_styles[2] = {"R4", "Q4"};
//};
//
//// float linear_basis_function(float, int);
//// float derivative_linear_basis_function(float, int, float);
//// float cubic_spline_basis_function(float, int);
//// float derivative_cubic_spline_basis_function(float, int, float);
//// float bernstein_quadratic_basis_function(float, int);
//// float derivative_bernstein_quadratic_basis_function(float, int, float);
//
//
//#endif
//#endif
