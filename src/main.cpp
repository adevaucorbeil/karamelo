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
#include <Kokkos_Core.hpp>

/*! Main program to drive MPM. */

int main(int argc, char **argv) {
  Kokkos::initialize();
  MPI_Init(&argc, &argv);                         /// Initialized MPI

  MPM(argc, argv, MPI_COMM_WORLD).input->file(); /// Create the MPM entity, read input file and execute commands

  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();                                 /// Finalize MPI
  Kokkos::finalize();
}
