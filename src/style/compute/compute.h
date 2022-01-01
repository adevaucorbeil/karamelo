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

#ifndef MPM_COMPUTE_H
#define MPM_COMPUTE_H

#include <pointers.h>
#include <vector>

class Compute : protected Pointers {
 public:
  string id;
  int igroup, groupbit;

  Compute(class MPM *, vector<string>);
  virtual ~Compute() {};
  virtual void init() = 0;
  virtual void setup() = 0;
  
  virtual void compute_value() = 0;
};

#endif
