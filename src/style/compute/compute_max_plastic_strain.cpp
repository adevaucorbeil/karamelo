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

#include <compute_max_plastic_strain.h>
#include <domain.h>
#include <error.h>
#include <group.h>
#include <input.h>
#include <math_special.h>
#include <output.h>
#include <universe.h>
#include <update.h>
#include <matrix.h>
#include <iostream>
#include <vector>

using namespace std;
using namespace MathSpecial;


ComputeMaxPlasticStrain::ComputeMaxPlasticStrain(MPM *mpm, vector<string> args)
    : Compute(mpm, args) {
  if (args.size() < 3) {
    error->all(FLERR, "Error: too few arguments for compute_max_plastic_strain: "
                      "requires at least 3 arguments. " +
                          to_string(args.size()) + " received\n");
  }

  if (group->pon[igroup] != "particles" &&
      group->pon[igroup] != "all") {
    error->all(
        FLERR,
        "compute_max_plastic_strain needs to be given a group of particles " +
            group->pon[igroup] + ", " + args[2] + " is a group of " +
            group->pon[igroup] + ".\n");
  }

  
  if (universe->me == 0)
    cout << "Creating new compute ComputeMaxPlasticStrain with ID: " << args[0] << endl;
  id = args[0];

  (*input->vars)[id + "_Epmax"] = Var(id + "_Epmax", 0);
  (*input->vars)[id + "_Tmax"] = Var(id + "_Tmax", 0);
}

void ComputeMaxPlasticStrain::compute_value(Solid &solid) {

  double Epmax(0.), Epmax_reduced(0.), Tmax(0.), Tmax_reduced(0.);

  Kokkos::View<double*> T = solid.T;
  Kokkos::View<double*> eff_plastic_strain = solid.eff_plastic_strain;
  Kokkos::View<int*> mask = solid.mask;

  int groupbit = this->groupbit;

  if (update->ntimestep == output->next ||
      update->ntimestep == update->nsteps)
    Kokkos::parallel_reduce("ComputeMaxPlasticStrain::compute_value", solid.np_local,
			    KOKKOS_LAMBDA(const int &ip, double &lTmax, double &lEpmax) {
			      if (mask[ip] & groupbit) {
				lTmax = lTmax > T[ip] ? lTmax : T[ip];
				lEpmax = lEpmax > eff_plastic_strain[ip] ? lEpmax : eff_plastic_strain[ip];
			      }
			    },Kokkos::Max<double>(Tmax), Kokkos::Max<double>(Epmax));

  // Reduce Epmax:
  MPI_Allreduce(&Epmax, &Epmax_reduced, 1, MPI_DOUBLE, MPI_MAX, universe->uworld);
  (*input->vars)[id + "_Epmax"] = Var(id + "_Epmax", Epmax_reduced);

  // Reduce Tmax:
  MPI_Allreduce(&Tmax, &Tmax_reduced, 1, MPI_DOUBLE, MPI_MAX, universe->uworld);
  (*input->vars)[id + "_Tmax"] = Var(id + "_Tmax", Tmax_reduced);
}
