/* -*- c++ -*- ----------------------------------------------------------*/

#ifndef MPM_LOG_H
#define MPM_LOG_H

#include "pointers.h"
#include <vector>

class Log : protected Pointers {
public:
  string style;

  Log(class MPM*, vector<string>);
  ~Log() {};

  void write();
};

#endif
