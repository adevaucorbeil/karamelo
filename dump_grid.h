/* -*- c++ -*- ----------------------------------------------------------*/

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
