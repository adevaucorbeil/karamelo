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

#ifdef COMMAND_CLASS

CommandStyle(internal_force,IntForce)

#else

#ifndef LMP_INTERNAL_FORCE_H
#define LMP_INTERNAL_FORCE_H

#include <pointers.h>

class IntForce : protected Pointers {
 public:
  IntForce(class MPM *);
  class Var command(vector<string>);

 private:
  int igroup, dir;
};

#endif
#endif
