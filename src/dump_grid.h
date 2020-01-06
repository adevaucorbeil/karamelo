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

DumpStyle(grid,DumpGrid)

#else

#ifndef MPM_DUMP_GRID_H
#define MPM_DUMP_GRID_H

#include "dump.h"

class DumpGrid : public Dump {
 public:
  DumpGrid(MPM *, vector<string>);
  ~DumpGrid();

  void write();
 protected:
  vector<string> known_var = {"x", "y", "z",
			      "vx", "vy", "vz",
			      "bx", "by", "bz",
			      "mass",
			      "ntypex", "ntypey", "ntypez"};
};

#endif
#endif
