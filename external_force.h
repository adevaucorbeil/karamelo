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

CommandStyle(external_force,ExtForce)

#else

#ifndef LMP_EXTERNAL_FORCE_H
#define LMP_EXTERNAL_FORCE_H

#include "pointers.h"

class ExtForce : protected Pointers {
 public:
  ExtForce(class MPM *);
  class Var command(vector<string>);

 private:
  int igroup, dir;
};

#endif
#endif
