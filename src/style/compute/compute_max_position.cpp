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

#include <compute_max_position.h>
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


ComputeMaxPosition::ComputeMaxPosition(MPM *mpm, vector<string> args)
    : Compute(mpm, args) {
  if (args.size() < 3) {
    error->all(FLERR, "Error: too few arguments for compute_max_position: "
                      "requires at least 3 arguments. " +
                          to_string(args.size()) + " received\n");
  }

  if (group->pon[igroup] != "particles" &&
      group->pon[igroup] != "all") {
    error->all(
        FLERR,
        "compute_max_position needs to be given a group of particles " +
            group->pon[igroup] + ", " + args[2] + " is a group of " +
            group->pon[igroup] + ".\n");
  }

  
  if (universe->me == 0)
    cout << "Creating new compute ComputeMaxPosition with ID: " << args[0] << endl;
  id = args[0];

  input->parsev(id + "_x", 0);
  input->parsev(id + "_y", 0);
  input->parsev(id + "_z", 0);

  (*input->vars)[id + "_x"] = Var(id + "_x", 0);
  (*input->vars)[id + "_y"] = Var(id + "_y", 0);
  (*input->vars)[id + "_z"] = Var(id + "_z", 0);
  t = update->ntimestep;
  Xmax[0] = Xmax[1] = Xmax[2] = 0;
}

void ComputeMaxPosition::compute_value(Solid &solid) {

  if (t != update->ntimestep) {
    t = update->ntimestep;
    Xmax[0] = Xmax[1] = Xmax[2] = 0;
  }

  float Xmax_reduced[3] = {0., 0., 0.};
  float Xmax_tmp[3] = {0., 0., 0.};
  Kokkos::View<Vector3d*> sx = solid.x;
  Kokkos::View<int*> mask = solid.mask;

  int groupbit = this->groupbit;

  if (update->ntimestep == output->next ||
      update->ntimestep == update->nsteps) {
    Kokkos::parallel_reduce("ComputeMaxPosition::compute_value", solid.np_local,
    KOKKOS_LAMBDA(const int &ip, float &lXmax0, float &lXmax1, float &lXmax2)
    {
      if (mask[ip] & groupbit) {
	const Vector3d &x = sx[ip];
	lXmax0 = lXmax0 > x[0] ? lXmax0 : x[0];
	lXmax1 = lXmax1 > x[1] ? lXmax1 : x[1];
	lXmax2 = lXmax2 > x[2] ? lXmax2 : x[2];
      }
    },Kokkos::Max<float>(Xmax_tmp[0]),
      Kokkos::Max<float>(Xmax_tmp[1]),
      Kokkos::Max<float>(Xmax_tmp[2]));

    for (int i = 0; i < domain->dimension; i++)
      Xmax[i] = max(Xmax[i], Xmax_tmp[i]);
  }

  // Reduce:
  MPI_Allreduce(Xmax, Xmax_reduced, 3, MPI_FLOAT, MPI_MAX, universe->uworld);
  input->parsev(id + "_x", Xmax_reduced[0]);
  input->parsev(id + "_y", Xmax_reduced[1]);
  input->parsev(id + "_z", Xmax_reduced[2]);
  (*input->vars)[id + "_x"] = Var(id + "_x", Xmax_reduced[0]);
  (*input->vars)[id + "_y"] = Var(id + "_y", Xmax_reduced[1]);
  (*input->vars)[id + "_z"] = Var(id + "_z", Xmax_reduced[2]);
}
