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
   memory(ptr->memory), 
   input(ptr->input), 
   output(ptr->output),
   infile(&ptr->infile),
   logfile(&ptr->logfile),
   domain(ptr->domain),
   material(ptr->material),
   variables(&ptr->variables),
   update(ptr->update),
   modify(ptr->modify),
   group(ptr->group) {}
  virtual ~Pointers() {}
 protected:
  MPM *mpm;
  Memory *&memory;
  Input *&input;
  Output *&output;
  
  Domain *&domain;
  Material *&material;
  Update *&update;
  Modify *&modify;
  Group *&group;

  filebuf *infile;
  filebuf *logfile;

  map<string, double> *variables;
};

  
#endif
