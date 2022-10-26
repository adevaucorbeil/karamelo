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

#include <compute_average_temperature.h>
#include <domain.h>
#include <error.h>
#include <group.h>
#include <input.h>
#include <math_special.h>
#include <output.h>
#include <universe.h>
#include <update.h>
#include <iostream>
#include <string>
#include <vector>

using namespace std;
using namespace MathSpecial;

ComputeAverageTemperature::ComputeAverageTemperature(MPM *mpm, vector<string> args)
    : Compute(mpm, args) {

  if (args.size() < 3) {
    error->all(FLERR, "Error: too few arguments for compute_average_temperature: "
                      "requires at least 3 arguments. " +
                          to_string(args.size()) + " received.\n");
  }

  if (group->pon[igroup] != "particles" &&
      group->pon[igroup] != "all") {
    error->all(FLERR,
               "compute_average_temperature needs to be given a group of particles" +
                   group->pon[igroup] + ", " + args[2] + " is a group of " +
                   group->pon[igroup] + ".\n");
  }

  
  if (universe->me == 0) 
    cout << "Creating new compute ComputeAverageTemperature with ID: " << args[0] << endl;
  id = args[0];

  input->parsev(id, 0);
}

ComputeAverageTemperature::~ComputeAverageTemperature() {}

void ComputeAverageTemperature::compute_value(Solid &solid) {

  float T, T_reduced;

  Kokkos::View<float*> sT = solid.T;
  Kokkos::View<int*> mask = solid.mask;

  int groupbit = this->groupbit;

  if (update->ntimestep == output->next ||
      update->ntimestep == update->nsteps)

  Kokkos::parallel_reduce("ComputeAverageTemperature::compute_value", solid.np_local,
    			  KOKKOS_LAMBDA(const int &ip, float &lT) {
			      if (mask[ip] & groupbit)
				lT += sT[ip];
			  },T);


  // Reduce T:
  MPI_Allreduce(&T, &T_reduced, 1, MPI_FLOAT, MPI_SUM, universe->uworld);
  input->parsev(id, T_reduced/group->n_tot(igroup));
}
