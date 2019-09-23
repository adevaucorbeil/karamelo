/* -*- c++ -*- ----------------------------------------------------------*/

#ifdef DUMP_CLASS

DumpStyle(grid/gz,DumpGridGz)

#else

#ifndef MPM_DUMP_GRID_GZ_H
#define MPM_DUMP_GRID_GZ_H

#include "dump.h"

class DumpGridGz : public Dump {
 public:
  DumpGridGz(MPM *, vector<string>);
  ~DumpGridGz();

  void write();
  //protected:
};

#endif
#endif
