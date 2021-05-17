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

#include "compute_average_velocity.h"
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
using namespace KARAMELO_NS;

ComputeAverageVelocity::ComputeAverageVelocity(MPM *mpm, vector<string> args)
    : Compute(mpm, args) {

  if (args.size() < 3) {
    error->all(FLERR, "Error: too few arguments for compute_average_velocity: "
                      "requires at least 3 arguments. " +
                          to_string(args.size()) + " received.\n");
  }

  if (group->pon[igroup].compare("particles") != 0 &&
      group->pon[igroup].compare("all") != 0) {
    error->all(FLERR,
               "compute_average_velocity needs to be given a group of particles" +
                   group->pon[igroup] + ", " + args[2] + " is a group of " +
                   group->pon[igroup] + ".\n");
  }

  cout << "Creating new compute ComputeAverageVelocity with ID: " << args[0] << endl;
  id = args[0];

  (*input->vars)[id + "_x"]=Var(id + "_x", 0);
  (*input->vars)[id + "_y"]=Var(id + "_y", 0);
  (*input->vars)[id + "_z"]=Var(id + "_z", 0);
}

ComputeAverageVelocity::~ComputeAverageVelocity() {}

void ComputeAverageVelocity::init() {}

void ComputeAverageVelocity::setup() {}

void ComputeAverageVelocity::compute_value() {
  if (update->ntimestep != output->next && update->ntimestep != update->nsteps)
    return;
  // cout << "In ComputeAverageVelocity::post_particles_to_grid()\n";

  // Go through all the nodes in the group and set b to the right value:
  double ux, uy, uz;

  int solid = group->solid[igroup];

  Solid *s;

  double vx, vy, vz, vx_reduced, vy_reduced, vz_reduced;
  int n = 0, n_reduced = 0;

  vx = vy = vz;
  vx_reduced = vz_reduced = vy_reduced = 0;

  if (solid == -1) {
    for (int isolid = 0; isolid < domain->solids.size(); isolid++) {
      s = domain->solids[isolid];

      for (int in = 0; in < s->np_local; in++) {
        if (s->mask[in] & groupbit) {
          vx += s->v[in](0);
          vy += s->v[in](1);
          vz += s->v[in](2);
	  n++;
        }
      }
    }
  } else {
    s = domain->solids[solid];

    for (int in = 0; in < s->np_local; in++) {
      if (s->mask[in] & groupbit) {
	vx += s->v[in](0);
	vy += s->v[in](1);
	vz += s->v[in](2);
	n++;
      }
    }
  }

  // Reduce Ek:
  MPI_Allreduce(&vx, &vx_reduced, 1, MPI_DOUBLE, MPI_SUM, universe->uworld);
  MPI_Allreduce(&vy, &vy_reduced, 1, MPI_DOUBLE, MPI_SUM, universe->uworld);
  MPI_Allreduce(&vz, &vz_reduced, 1, MPI_DOUBLE, MPI_SUM, universe->uworld);
  MPI_Allreduce(&n, &n_reduced, 1, MPI_INT, MPI_SUM, universe->uworld);

  (*input->vars)[id + "_x"]=Var(id + "_x", vx_reduced/n);
  (*input->vars)[id + "_y"]=Var(id + "_y", vy_reduced/n);
  (*input->vars)[id + "_z"]=Var(id + "_z", vz_reduced/n);
}
