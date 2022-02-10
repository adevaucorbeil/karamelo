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

ComputeStyle(average_stress,ComputeAverageStress)

#else

#ifndef MPM_COMPUTE_AVERAGE_STRESS_H
#define MPM_COMPUTE_AVERAGE_STRESS_H

#include <compute.h>
#include <var.h>
#include <vector>

//namespace KARAMELO_NS {

class ComputeAverageStress : public Compute {
 public:
  ComputeAverageStress(class MPM *, vector<string>);
  ~ComputeAverageStress();
  void init();
  void setup();

  void compute_value();

private:
  // class Var xvalue, yvalue, zvalue;    // Set force in x, y, and z directions.
  // bool xset, yset, zset;               // Does the compute set the x, y, and z forces of the group?
};
//}
#endif
#endif

