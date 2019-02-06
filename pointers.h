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
   output(ptr->output),
   infile(&ptr->infile),
   logfile(&ptr->logfile),
   domain(ptr->domain),
   material(ptr->material),
   variables(&ptr->variables),
   update(ptr->update),
   modify(ptr->modify) {}
  virtual ~Pointers() {}
 protected:
  MPM *mpm;
  Input *&input;
  Output *&output;
  
  Domain *&domain;
  Material *&material;
  Update *&update;
  Modify *&modify;

  filebuf *infile;
  filebuf *logfile;

  map<string, double> *variables;
};

  
#endif
