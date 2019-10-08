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

#ifndef MPM_SCHEME_H
#define MPM_SCHEME_H

#include "pointers.h"

class Scheme : protected Pointers {
 public:
  Scheme(class MPM *);
  virtual ~Scheme();
  virtual void setup() = 0;
  virtual void run(class Var) = 0;
};

#endif
