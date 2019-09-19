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
};

#endif
#endif
