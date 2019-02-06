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
  void options(vector<string> *, vector<string>::iterator);

  // implemented by each dump

 protected:
  string filename;
};

#endif
