/* -*- c++ -*- ----------------------------------------------------------*/


#ifdef COMMAND_CLASS

CommandStyle(run,Run)

#else

#ifndef MPM_RUN_H
#define MPM_RUN_H

#include "pointers.h"

class Run : protected Pointers {
 public:
  Run(class MPM *);
  void command(vector<string>);
};

#endif
#endif
