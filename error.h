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

#ifndef MPM_ERROR_H
#define MPM_ERROR_H

#include "pointers.h"

class Error : protected Pointers {
 public:
  Error(class MPM *);

  void all(const char *, int, const string);
  void one(const char *, int, const string);
  void done(int = 0); // 1 would be fully backwards compatible
};

#endif
