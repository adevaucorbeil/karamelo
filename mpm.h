/* -*- c++ -*- ----------------------------------------------------------*/

#ifndef MPM_MPM_H
#define MPM_MPM_H

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <string>
#include <fstream>

using namespace std;

class MPM {

 public:

  class Input *input;            // input script processing
  class Solid *solid;            // solid list
  //class *settings; // settings
  filebuf infile;                // infile
  filebuf logfile;               // logfile

  MPM(int, char **);
  ~MPM();
};

#endif
