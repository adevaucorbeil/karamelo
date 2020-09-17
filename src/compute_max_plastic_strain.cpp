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

#include "compute_max_plastic_strain.h"
#include "domain.h"
#include "error.h"
#include "group.h"
#include "input.h"
#include "math_special.h"
#include "output.h"
#include "universe.h"
#include "update.h"
#include <Eigen/Eigen>
#include <iostream>
#include <vector>

using namespace std;
using namespace MathSpecial;
using namespace Eigen;

ComputeMaxPlasticStrain::ComputeMaxPlasticStrain(MPM *mpm, vector<string> args)
    : Compute(mpm, args) {
  if (args.size() < 3) {
    error->all(FLERR, "Error: too few arguments for compute_max_plastic_strain: "
                      "requires at least 3 arguments. " +
                          to_string(args.size()) + " received\n");
  }

  if (group->pon[igroup].compare("particles") != 0 &&
      group->pon[igroup].compare("all") != 0) {
    error->all(
        FLERR,
        "compute_max_plastic_strain needs to be given a group of particles " +
            group->pon[igroup] + ", " + args[2] + " is a group of " +
            group->pon[igroup] + ".\n");
  }

  cout << "Creating new compute ComputeMaxPlasticStrain with ID: " << args[0] << endl;
  id = args[0];

  (*input->vars)[id + "_Epmax"] = Var(id + "_Epmax", 0);
  (*input->vars)[id + "_Tmax"] = Var(id + "_Tmax", 0);
}

ComputeMaxPlasticStrain::~ComputeMaxPlasticStrain() {}

void ComputeMaxPlasticStrain::init() {}

void ComputeMaxPlasticStrain::setup() {}

void ComputeMaxPlasticStrain::compute_value() {
  if (update->ntimestep != output->next && update->ntimestep != update->nsteps)
    return;
  // cout << "In ComputeMaxPlasticStrain::post_particles_to_grid()\n";

  // Go through all the nodes in the group and set b to the right value:
  double ux, uy, uz;

  int solid = group->solid[igroup];

  Solid *s;

  double Epmax(0.), Epmax_reduced(0.), Tmax(0.), Tmax_reduced(0.);

  if (solid == -1) {
    for (int isolid = 0; isolid < domain->solids.size(); isolid++) {
      s = domain->solids[isolid];

      for (int in = 0; in < s->np_local; in++) {
        if (s->mask[in] & groupbit) {
          Epmax = max(Epmax, s->eff_plastic_strain[in]);
          Tmax = max(Tmax, s->T[in]);
        }
      }
    }
  } else {
    s = domain->solids[solid];

    for (int in = 0; in < s->np_local; in++) {
      if (s->mask[in] & groupbit) {
	Epmax = max(Epmax, s->eff_plastic_strain[in]);
	Tmax = max(Tmax, s->T[in]);
      }
    }

  }

  // Reduce Epmax:
  MPI_Allreduce(&Epmax, &Epmax_reduced, 1, MPI_DOUBLE, MPI_MAX, universe->uworld);
  (*input->vars)[id + "_Epmax"] = Var(id + "_Epmax", Epmax_reduced);

  // Reduce Tmax:
  MPI_Allreduce(&Tmax, &Tmax_reduced, 1, MPI_DOUBLE, MPI_MAX, universe->uworld);
  (*input->vars)[id + "_Tmax"] = Var(id + "_Tmax", Tmax_reduced);
}
