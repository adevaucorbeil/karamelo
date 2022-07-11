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

#include <compute_average_velocity.h>
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
#include <string>
#include <vector>

using namespace std;
using namespace MathSpecial;

//using namespace KARAMELO_NS;

ComputeAverageVelocity::ComputeAverageVelocity(MPM *mpm, vector<string> args)
    : Compute(mpm, args) {

  if (args.size() < 3) {
    error->all(FLERR, "Error: too few arguments for compute_average_velocity: "
                      "requires at least 3 arguments. " +
                          to_string(args.size()) + " received.\n");
  }

  if (group->pon[igroup] != "particles" &&
      group->pon[igroup] != "all") {
    error->all(FLERR,
               "compute_average_velocity needs to be given a group of particles" +
                   group->pon[igroup] + ", " + args[2] + " is a group of " +
                   group->pon[igroup] + ".\n");
  }

  if (universe->me == 0)
    cout << "Creating new compute ComputeAverageVelocity with ID: " << args[0] << endl;
  id = args[0];

  t = update->ntimestep;
  input->parsev(id + "_x", 0);
  input->parsev(id + "_y", 0);
  input->parsev(id + "_z", 0);

  v_average[0] = v_average[1] = v_average[2] = 0;
}

ComputeAverageVelocity::~ComputeAverageVelocity() {}

void ComputeAverageVelocity::compute_value(Solid &solid) {
  if (t != update->ntimestep) {
    t = update->ntimestep;
    v_average[0] = v_average[1] = v_average[2] = 0;
  }

  float vx, vy, vz;

  vx = vy = vz = 0;

  Kokkos::View<Vector3d*> sv = solid.v;
  Kokkos::View<int*> mask = solid.mask;

  int groupbit = this->groupbit;

  if (update->ntimestep == output->next ||
      update->ntimestep == update->nsteps) {
    Kokkos::parallel_reduce("ComputeAverageVelocity::compute_value", solid.np_local,
			    KOKKOS_LAMBDA(const int &ip,
					  float &lvx, float &lvy, float &lvz) {
			      if(mask[ip] & groupbit) {
				lvx = sv[ip](0);
				lvy = sv[ip](1);
				lvz = sv[ip](2);
			      }
			    },vx, vy, vz);

    v_average[0] += vx;
    v_average[1] += vy;
    v_average[2] += vz;
  }


  // Reduce velocity:
  float v_reduced[3] = {0, 0, 0};
  MPI_Allreduce(v_average, v_reduced, 3, MPI_FLOAT, MPI_SUM, universe->uworld);

  input->parsev(id + "_x", v_reduced[0]/group->n_tot(igroup));
  input->parsev(id + "_y", v_reduced[1]/group->n_tot(igroup));
  input->parsev(id + "_z", v_reduced[2]/group->n_tot(igroup));
}
