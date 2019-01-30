/* -*- c++ -*- ----------------------------------------------------------*/

#ifndef MPM_MPM_H
#define MPM_MPM_H

#include <stdio.h>
#include <stdlib.h>

class MPM {

 public:

  class Input *input;            // input script processing
  class Solid *solid;            // solid list
  //class *settings; // settings
  FILE *infile;                  // infile
  FILE *logfile;                 // logfile

  MPM(int, char **);
  ~MPM();
};

#endif
