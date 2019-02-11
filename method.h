/* -*- c++ -*- ----------------------------------------------------------*/

#ifndef MPM_METHOD_H
#define MPM_METHOD_H

#include "pointers.h"

class Method : protected Pointers {
 public:
  Method(class MPM *);
  virtual ~Method();

  virtual void setup() = 0;
  virtual void compute_grid_weight_functions_and_gradients() = 0;
  virtual void particles_to_grid() = 0;
  virtual void update_grid_state() = 0;
  virtual void grid_to_points() = 0;
  virtual void advance_particles() = 0;
  virtual void velocities_to_grid() = 0;
  virtual void compute_rate_deformation_gradient() = 0;
  virtual void update_deformation_gradient() = 0;
  virtual void update_stress() = 0;
};

#endif
