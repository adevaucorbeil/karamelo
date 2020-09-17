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

#include "compute_kinetic_energy.h"
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
#include <string>
#include <vector>

using namespace std;
using namespace MathSpecial;
using namespace Eigen;

ComputeKineticEnergy::ComputeKineticEnergy(MPM *mpm, vector<string> args)
    : Compute(mpm, args) {

  if (args.size() < 3) {
    error->all(FLERR, "Error: too few arguments for compute_kinetic_energy: "
                      "requires at least 3 arguments. " +
                          to_string(args.size()) + " received.\n");
  }

  if (group->pon[igroup].compare("particles") != 0 &&
      group->pon[igroup].compare("all") != 0) {
    error->all(FLERR,
               "compute_kinetic_energy needs to be given a group of particles" +
                   group->pon[igroup] + ", " + args[2] + " is a group of " +
                   group->pon[igroup] + ".\n");
  }

  cout << "Creating new compute ComputeKineticEnergy with ID: " << args[0] << endl;
  id = args[0];

  (*input->vars)[id]=Var(id, 0);
}

ComputeKineticEnergy::~ComputeKineticEnergy() {}

void ComputeKineticEnergy::init() {}

void ComputeKineticEnergy::setup() {}

void ComputeKineticEnergy::compute_value() {
  if (update->ntimestep != output->next && update->ntimestep != update->nsteps)
    return;
  // cout << "In ComputeKineticEnergy::post_particles_to_grid()\n";

  // Go through all the nodes in the group and set b to the right value:
  double ux, uy, uz;

  int solid = group->solid[igroup];

  Solid *s;

  double Ek, Ek_reduced;

  Ek = 0;
  Ek_reduced = 0;

  if (solid == -1) {
    for (int isolid = 0; isolid < domain->solids.size(); isolid++) {
      s = domain->solids[isolid];

      for (int in = 0; in < s->np_local; in++) {
        if (s->mask[in] & groupbit) {
          Ek += 0.5 * s->mass[in] * square(s->v[in].norm());
        }
      }
    }
  } else {
    s = domain->solids[solid];

    for (int in = 0; in < s->np_local; in++) {
      if (s->mask[in] & groupbit) {
        Ek += 0.5 * s->mass[in] * square(s->v[in].norm());
      }
    }
  }

  // Reduce Ek:
  MPI_Allreduce(&Ek, &Ek_reduced, 1, MPI_DOUBLE, MPI_SUM, universe->uworld);
  (*input->vars)[id] = Var(id, Ek_reduced);
}
