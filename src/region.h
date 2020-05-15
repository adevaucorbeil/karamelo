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

/* \page region
 */

#ifndef MPM_REGION_H
#define MPM_REGION_H

#include "pointers.h"
#include <vector>

class Region : protected Pointers {
 public:
  string id;
  int interior;                     // 1 for interior, 0 for exterior

  Region(class MPM *, vector<string>);
  virtual ~Region();
  virtual void init();
  void options(vector<string> *, vector<string>::iterator);

  // called by other classes to check point versus region

  int match(double, double, double);

  // implemented by each region

  virtual vector<double> limits() {return vector<double>();};
  virtual int inside(double, double, double) = 0;
  //protected:
};

#endif
