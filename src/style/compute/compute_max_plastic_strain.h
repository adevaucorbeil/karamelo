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

ComputeStyle(max_plastic_strain, ComputeMaxPlasticStrain)

#else

#ifndef MPM_COMPUTE_MAX_PLASTIC_STRAIN_H
#define MPM_COMPUTE_MAX_PLASTIC_STRAIN_H

#include <compute.h>
#include <vector>

class ComputeMaxPlasticStrain : public Compute {
public:
  ComputeMaxPlasticStrain(class MPM *, vector<string>);

  void compute_value(Solid &solid);

private:
  int t;
  float Epmax, Tmax;
};

#endif
#endif
