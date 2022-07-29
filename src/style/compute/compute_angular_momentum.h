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

ComputeStyle(angular_momentum,ComputeAngularMomentum)

#else

#ifndef MPM_COMPUTE_ANGULAR_MOMENTUM_H
#define MPM_COMPUTE_ANGULAR_MOMENTUM_H

#include <compute.h>
#include <matrix.h>

class ComputeAngularMomentum : public Compute {
 public:
  ComputeAngularMomentum(class MPM *, vector<string>);
  ~ComputeAngularMomentum();

  void compute_value(Solid &solid);

private:
  string usage = "Usage: compute(computer-ID, angular_momentum, group, x0, y0, z0)\n";
  int Nargs = 6;
  int t;
  Vector3d J, x0;
};

#endif
#endif

