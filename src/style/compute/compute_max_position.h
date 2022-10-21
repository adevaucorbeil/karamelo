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

ComputeStyle(max_position, ComputeMaxPosition)

#else

#ifndef MPM_COMPUTE_MAX_POSITION_H
#define MPM_COMPUTE_MAX_POSITION_H

#include <compute.h>
#include <vector>

class ComputeMaxPosition : public Compute {
public:
  ComputeMaxPosition(class MPM *, vector<string>);

  void compute_value(Solid &solid);

private:
  int t;
  float Xmax[3];
};

#endif
#endif
