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

SchemeStyle(musl,MUSL)

#else

#ifndef LMP_MUSL_H
#define LMP_MUSL_H

#include <scheme.h>
#include <var.h>
#include <vector>

class MUSL : public Scheme {
 public:
  MUSL(class MPM *);
  ~MUSL() {}
  void setup();
  void run(class Var);
};

#endif
#endif
