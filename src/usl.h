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

#ifdef SCHEME_CLASS

SchemeStyle(usl,USL)

#else

#ifndef LMP_USL_H
#define LMP_USL_H

#include "scheme.h"
#include "var.h"
#include <vector>

class USL : public Scheme {
 public:
  USL(class MPM *);
  ~USL() {}
  void setup();
  void run(class Var);
};

#endif
#endif
