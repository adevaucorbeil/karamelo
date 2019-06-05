/* -*- c++ -*- ----------------------------------------------------------*/


#ifdef COMMAND_CLASS

CommandStyle(run_time,RunTime)

#else

#ifndef MPM_RUN_TIME_H
#define MPM_RUN_TIME_H

#include "pointers.h"

class RunTime : protected Pointers {
 public:
  RunTime(class MPM *);
  class Var command(vector<string>);
};

#endif
#endif
