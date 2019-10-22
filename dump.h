/* -*- c++ -*- ----------------------------------------------------------*/

#ifndef MPM_DUMP_H
#define MPM_DUMP_H

#include "pointers.h"
#include <vector>

class Dump : protected Pointers {
 public:
  string id;
  string style;
  //int group; //groups are not supported yet: default all
  

  Dump(class MPM *, vector<string>);
  virtual ~Dump();

  // implemented by each dump

  virtual void write() = 0;

 protected:
  string filename;
  vector<string> output_var;
};

#endif
