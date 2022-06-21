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

#ifdef COMPUTE_CLASS

ComputeStyle(average_velocity,ComputeAverageVelocity)

#else

#ifndef MPM_COMPUTE_AVERAGE_VELOCITY_H
#define MPM_COMPUTE_AVERAGE_VELOCITY_H

#include <compute.h>
#include <vector>

class ComputeAverageVelocity : public Compute {
 public:
  ComputeAverageVelocity(class MPM *, vector<string>);
  ~ComputeAverageVelocity();

  void compute_value(Solid &solid);

private:
  int t;
  double v_average[3];
};

#endif
#endif

