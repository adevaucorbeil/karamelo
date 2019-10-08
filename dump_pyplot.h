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

#ifdef DUMP_CLASS

DumpStyle(pyplot,DumpPyPlot)

#else

#ifndef MPM_DUMP_PYPLOT_H
#define MPM_DUMP_PYPLOT_H

#include "dump.h"

class DumpPyPlot : public Dump {
 public:
  DumpPyPlot(MPM *, vector<string>);
  ~DumpPyPlot();

  void write();
  //protected:
 private:
  int xsize, ysize;
};

#endif
#endif
