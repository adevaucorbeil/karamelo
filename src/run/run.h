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

CommandStyle(run,Run)

#else

#ifndef MPM_RUN_H
#define MPM_RUN_H

#include <pointers.h>

class Run : protected Pointers {
 public:
  Run(class MPM *);
  class Var command(vector<string>);
};

#endif
#endif
