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

#ifndef LMP_MATH_SPECIAL_H
#define LMP_MATH_SPECIAL_H

#include <math.h>


namespace MathSpecial {
  // x**2, use instead of pow(x,2.0)

  static inline double square(const double &x) { return x*x; }
}

#endif
