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

CommandStyle(run_while,RunWhile)

#else

#ifndef MPM_RUN_WHILE_H
#define MPM_RUN_WHILE_H

#include "pointers.h"

class RunWhile : protected Pointers {
 public:
  RunWhile(class MPM *);
  class Var command(vector<string>);
};

#endif
#endif
