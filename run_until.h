/* -*- c++ -*- ----------------------------------------------------------*/


#ifdef COMMAND_CLASS

CommandStyle(run_until,RunUntil)

#else

#ifndef MPM_RUN_UNTIL_H
#define MPM_RUN_UNTIL_H

#include "pointers.h"

class RunUntil : protected Pointers {
 public:
  RunUntil(class MPM *);
  class Var command(vector<string>);
};

#endif
#endif
