/* -*- c++ -*- ----------------------------------------------------------*/

#ifndef MPM_MPM_H
#define MPM_MPM_H

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <string>
#include <fstream>
#include <map>
#include "variable.h"

using namespace std;

class MPM {

 public:

  class Memory *memory;          // memory allocation functions
  class Input *input;            // input script processing
  class Output *output;          // thermo/dump/restart

  class Domain *domain;          // simulation box
  class Material *material;      // material
  class Update *update;          //
  class Modify *modify;          // fixes and computes
  class Group *group;            // groups of particles

  filebuf infile;                // infile
  filebuf logfile;               // logfile
  ofstream *wlogfile;

  map<string, class Variable> *variables; // global variables

  MPM(int, char **);
  ~MPM();
  void init();
};

#endif
