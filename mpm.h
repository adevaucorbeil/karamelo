/* -*- c++ -*- ----------------------------------------------------------*/

#ifndef MPM_MPM_H
#define MPM_MPM_H

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <string>
#include <fstream>
#include <map>

using namespace std;

class MPM {

 public:

  class Memory *memory;          // memory allocation functions
  class Universe *universe;      // universe of processors
  class Input *input;            // input script processing
  class Output *output;          // thermo/dump/restart

  class Domain *domain;          // simulation box
  class Material *material;      // material
  class Update *update;          //
  class Modify *modify;          // fixes and computes
  class Group *group;            // groups of particles

  MPI_Comm world;                // MPI communicator
  double initclock;              // wall clock at instantiation

  filebuf infile;                // infile
  filebuf logfile;               // logfile
  ofstream *wlogfile;

  MPM(int, char **);
  ~MPM();
  void init();
};

#endif
