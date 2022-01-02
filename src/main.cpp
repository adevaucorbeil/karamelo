/* ----------------------------------------------------------------------
 *
 *                    ***       Karamelo       ***
 *               Parallel Material Point Method Simulator
 * 
 * Copyright (2019) Alban de Vaucorbeil, alban.devaucorbeil@monash.edu
 * Materials Science and Engineering, Monash University
 * Clayton VIC 3800, Australia

 * This software is distributed under the GNU General Public License.
 *
 * ----------------------------------------------------------------------- */

#include <mpi.h>
#include <mpm.h>
#include <input.h>

/*! Main program to drive MPM. */

int main(int argc, char **argv) {
  MPI_Init(&argc,&argv);                        /// Initialized MPI

  MPM *mpm = new MPM(argc,argv,MPI_COMM_WORLD); /// Create the MPM entity
  mpm->input->file();                           /// Read input file and execute commands

  delete mpm;

  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();                               /// Finalize MPI
}
