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
   logfile(&ptr->logfile),
   domain(ptr->domain),
   variables(&ptr->variables) {}
  virtual ~Pointers() {}
 protected:
  MPM *mpm;
  Input *&input;

  
  Domain *&domain;
  Solid *&solid;

  filebuf *infile;
  filebuf *logfile;

  map<string, double> *variables;
};

  
#endif
