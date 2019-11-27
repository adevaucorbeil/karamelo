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

#ifndef MPM_DUMP_H
#define MPM_DUMP_H

#include "pointers.h"
#include <vector>

class Dump : protected Pointers {
 public:
  string id;
  string style;
  //int group; //groups are not supported yet: default all
  

  Dump(class MPM *, vector<string>);
  virtual ~Dump();

  // implemented by each dump

  virtual void write() = 0;

 protected:
  string filename;
  vector<string> output_var;
};

#endif
