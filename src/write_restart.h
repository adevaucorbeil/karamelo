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

CommandStyle(write_restart,WriteRestart)

#else

#ifndef MPM_WRITE_RESTART_H
#define MPM_WRITE_RESTART_H

#include "pointers.h"

class WriteRestart : protected Pointers {
 public:
  WriteRestart(class MPM *);
  class Var command(vector<string>);
  void write();

protected:
  string filename;
};

#endif
#endif
