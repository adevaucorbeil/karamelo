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

ComputeStyle(average_temperature,ComputeAverageTemperature)

#else

#ifndef MPM_COMPUTE_AVERAGE_TEMPERATURE_H
#define MPM_COMPUTE_AVERAGE_TEMPERATURE_H

#include <compute.h>
#include <vector>

class ComputeAverageTemperature : public Compute {
 public:
  ComputeAverageTemperature(class MPM *, vector<string>);
  ~ComputeAverageTemperature();

  void compute_value(Solid &solid);
};

#endif
#endif

