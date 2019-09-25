#include <stdio.h>
#include <mpi.h>
#include "mpm.h"
#include "input.h"

/* ----------------------------------------------------------------------
                        main program to drive MPM 
------------------------------------------------------------------------- */

int main(int argc, char **argv) {

  MPI_Init(&argc,&argv);

  MPM *mpm = new MPM(argc,argv,MPI_COMM_WORLD);
  mpm->input->file();

  delete mpm;

  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();
}
