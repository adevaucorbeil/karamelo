/* -*- c++ -*- ----------------------------------------------------------*/

#ifdef METHOD_CLASS

MethodStyle(tlmpm,TLMPM)

#else

#ifndef LMP_TLMPM_H
#define LMP_TLMPM_H

#include "method.h"
#include <vector>

class TLMPM : public Method {
 public:

  map<int, int> **neigh_pn; // List of the nodes neighbouring a given particle for each solid
  map<int, int> **neigh_np; // List of the particles neighbouring a given node for each solid

  TLMPM(class MPM *, vector<string>);
  ~TLMPM();

  void setup();
  void compute_grid_weight_functions_and_gradients();
  void particles_to_grid();
  void update_grid_state();
  void grid_to_points();
  void advance_particles();
  void velocities_to_grid();
  void compute_rate_deformation_gradient();
  void update_deformation_gradient();
  void update_stress();

 protected:
};

#endif
#endif
