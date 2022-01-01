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

CommandStyle(read_restart,ReadRestart)

#else

#ifndef MPM_READ_RESTART_H
#define MPM_READ_RESTART_H

#include <pointers.h>
#include <fstream>
#include <vector>

class ReadRestart : protected Pointers {
 public:
  ReadRestart(class MPM *);
  ~ReadRestart();
  class Var command(vector<string>);

private:
  string filename;
  size_t pos_asterisk;
  ifstream *ifr;
  void header();

  int read_int();
  string read_string();
};

#endif
#endif
