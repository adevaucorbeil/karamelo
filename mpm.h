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

  class Input *input;            // input script processing

  class Domain *domain;          // simulation box
  class Solid *solid;            // solid list
  //class *settings; // settings

  filebuf infile;                // infile
  filebuf logfile;               // logfile

  map<string, double> variables; // global variables

  MPM(int, char **);
  ~MPM();
};

#endif
