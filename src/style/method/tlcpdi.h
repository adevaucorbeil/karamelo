///* -*- c++ -*- ----------------------------------------------------------*/
//
//#ifdef METHOD_CLASS
//
//MethodStyle(tlcpdi,TLCPDI)
//
//#else
//
//#ifndef LMP_TLCPDI_H
//#define LMP_TLCPDI_H
//
//#include <method.h>
//#include <vector>
//#include <matrix.h>
//
//
//class TLCPDI : public Method {
// public:
//
//  TLCPDI(class MPM *);
//  ~TLCPDI();
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
//  void update_grid_positions();
//  void compute_rate_deformation_gradient(bool);
//  void update_deformation_gradient();
//  void update_stress(bool);
//  void adjust_dt();
//  void reset();
//  void exchange_particles() {};
//
//  bool update_wf;
//
//  private:
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
