/* -*- c++ -*- ----------------------------------------------------------*/


#ifdef COMMAND_CLASS

CommandStyle(run_while,RunWhile)

#else

#ifndef MPM_RUN_WHILE_H
#define MPM_RUN_WHILE_H

#include "pointers.h"

class RunWhile : protected Pointers {
 public:
  RunWhile(class MPM *);
  class Var command(vector<string>);
};

#endif
#endif
