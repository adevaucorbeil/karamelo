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

ComputeStyle(kinetic_energy,ComputeKineticEnergy)

#else

#ifndef MPM_COMPUTE_KINETIC_ENERGY_H
#define MPM_COMPUTE_KINETIC_ENERGY_H

#include <compute.h>
#include <vector>

class ComputeKineticEnergy : public Compute {
 public:
  ComputeKineticEnergy(class MPM *, vector<string>);
  ~ComputeKineticEnergy();

  void compute_value(Solid &solid);

private:
  int t;
  double Ek;
};

#endif
#endif

