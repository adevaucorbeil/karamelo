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

#include <pointers.h>
#include <fstream>
#include <vector>

class WriteRestart : protected Pointers {
 public:
  WriteRestart(class MPM *);
  ~WriteRestart();
  class Var command(vector<string>);
  void write();

private:
  string filename;
  size_t pos_asterisk;
  ofstream *of;
  void header();

  template <typename T> void write_variable(int, T);
  void write_string(int flag, string value);
  //void write_int(int flag, int value);
};

#endif
#endif
