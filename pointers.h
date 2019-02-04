/* -*- c++ -*- ----------------------------------------------------------*/

#ifndef MPM_POINTERS_H
#define MPM_POINTERS_H

#include "mpm.h"
#include "mpmtype.h"

using namespace std;

class Pointers {
 public:
 Pointers(MPM *ptr) :
   mpm(ptr),
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

  filebuf *infile;
  filebuf *logfile;

  map<string, double> *variables;
};

  
#endif
