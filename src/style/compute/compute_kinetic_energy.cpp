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

#include <compute_kinetic_energy.h>
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


ComputeKineticEnergy::ComputeKineticEnergy(MPM *mpm, vector<string> args)
    : Compute(mpm, args) {

  if (args.size() < 3) {
    error->all(FLERR, "Error: too few arguments for compute_kinetic_energy: "
                      "requires at least 3 arguments. " +
                          to_string(args.size()) + " received.\n");
  }

  if (group->pon[igroup] != "particles" &&
      group->pon[igroup] != "all") {
    error->all(FLERR,
               "compute_kinetic_energy needs to be given a group of particles" +
                   group->pon[igroup] + ", " + args[2] + " is a group of " +
                   group->pon[igroup] + ".\n");
  }

  
  if (universe->me == 0) 
    cout << "Creating new compute ComputeKineticEnergy with ID: " << args[0] << endl;
  id = args[0];

  (*input->vars)[id]=Var(id, 0);
}

ComputeKineticEnergy::~ComputeKineticEnergy() {}

void ComputeKineticEnergy::compute_value(Solid &solid) {
  double Ek, Ek_reduced;

  Ek = 0;
  Ek_reduced = 0;

  Kokkos::View<Vector3d*> v = solid.v;
  Kokkos::View<double*> mass = solid.mass;
  Kokkos::View<int*> mask = solid.mask;

  int groupbit = this->groupbit;
  int nsteps = update->nsteps;
  bigint next = output->next;
  double ntimestep = update->ntimestep;

  Kokkos::parallel_reduce("ComputeAverageStress::compute_value", solid.np_local,
			  KOKKOS_LAMBDA(const int &ip, double &lEk) {
			 if ((ntimestep != next &&
			      ntimestep != nsteps) ||
			     !(mask[ip] & groupbit))
			   return;

			 lEk += 0.5 * mass[ip] * square(v[ip].norm());
			  },Ek);

  // Reduce Ek:
  MPI_Allreduce(&Ek, &Ek_reduced, 1, MPI_DOUBLE, MPI_SUM, universe->uworld);
  (*input->vars)[id] = Var(id, Ek_reduced);
}
