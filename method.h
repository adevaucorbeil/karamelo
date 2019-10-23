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

#ifndef MPM_METHOD_H
#define MPM_METHOD_H

#include "pointers.h"
#include <vector>

class Method : protected Pointers {
 public:
  Method(class MPM *);
  virtual ~Method();

  virtual void setup(vector<string>) = 0;
  virtual void compute_grid_weight_functions_and_gradients() = 0;
  virtual void particles_to_grid() = 0;
  virtual void update_grid_state() = 0;
  virtual void grid_to_points() = 0;
  virtual void advance_particles() = 0;
  virtual void velocities_to_grid() = 0;
  virtual void update_grid_positions() = 0;
  virtual void compute_rate_deformation_gradient() = 0;
  virtual void update_deformation_gradient() = 0;
  virtual void update_stress() = 0;
  virtual void adjust_dt() = 0;
  virtual void reset() = 0;
  virtual void exchange_particles() = 0;

  bool is_TL;         // true: the method is total Lagrangian; false: it is updated Lagrangian
  bool is_CPDI;       // true if the method is a CPDI-like
};

#endif
