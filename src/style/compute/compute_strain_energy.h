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

ComputeStyle(strain_energy,ComputeStrainEnergy)

#else

#ifndef MPM_COMPUTE_STRAIN_ENERGY_H
#define MPM_COMPUTE_STRAIN_ENERGY_H

#include <compute.h>
#include <vector>

class ComputeStrainEnergy : public Compute {
 public:
  ComputeStrainEnergy(class MPM *, vector<string>);
  
  void compute_value(Solid &solid);

private:
  int t;
  float Es;
};

#endif
#endif

