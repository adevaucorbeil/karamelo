/* -*- c++ -*- ----------------------------------------------------------*/

#ifndef MPM_POINTERS_H
#define MPM_POINTERS_H

#include "mpm.h"

using namespace std;

class Pointers {
 public:
 Pointers(MPM *ptr) :
   mpm(ptr),
   solid(ptr->solid),
   input(ptr->input), 
   infile(&ptr->infile),
   logfile(&ptr->logfile) {}
  virtual ~Pointers() {}
 protected:
  MPM *mpm;
  Solid *&solid;
  Input *&input;

  filebuf *infile;
  filebuf *logfile;
};

  
#endif
