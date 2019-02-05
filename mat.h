/* -*- c++ -*- ----------------------------------------------------------*/

#ifndef MPM_MAT_H
#define MPM_MAT_H

#include "pointers.h"
#include "eos.h"
#include <vector>

class Mat : protected Pointers {
 public:
  string id;
  class EOS * eos;

  Mat(class MPM *, vector<string>);
  virtual ~Mat();
  virtual void init();
  void options(vector<string> *, vector<string>::iterator);

  // implemented by each mat
  //protected:
};

#endif
