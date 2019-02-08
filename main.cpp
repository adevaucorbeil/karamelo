#include <stdio.h>
#include "mpm.h"
#include "input.h"

/* ----------------------------------------------------------------------
                        main program to drive MPM 
------------------------------------------------------------------------- */

int main(int argc, char **argv) {

  MPM *mpm = new MPM(argc,argv);
  mpm->input->file();

  delete mpm;
}
