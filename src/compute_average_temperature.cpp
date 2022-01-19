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

#include "compute_average_temperature.h"
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

ComputeAverageTemperature::ComputeAverageTemperature(MPM *mpm, vector<string> args)
    : Compute(mpm, args) {

  if (args.size() < 3) {
    error->all(FLERR, "Error: too few arguments for compute_average_temperature: "
                      "requires at least 3 arguments. " +
                          to_string(args.size()) + " received.\n");
  }

  if (group->pon[igroup].compare("particles") != 0 &&
      group->pon[igroup].compare("all") != 0) {
    error->all(FLERR,
               "compute_average_temperature needs to be given a group of particles" +
                   group->pon[igroup] + ", " + args[2] + " is a group of " +
                   group->pon[igroup] + ".\n");
  }

  
  if (universe->me == 0) 
    cout << "Creating new compute ComputeAverageTemperature with ID: " << args[0] << endl;
  id = args[0];

  (*input->vars)[id]=Var(id, 0);
}

ComputeAverageTemperature::~ComputeAverageTemperature() {}

void ComputeAverageTemperature::init() {}

void ComputeAverageTemperature::setup() {}

void ComputeAverageTemperature::compute_value() {
  if (update->ntimestep != output->next && update->ntimestep != update->nsteps)
    return;

  // Go through all the nodes in the group and set b to the right value:

  int solid = group->solid[igroup];

  Solid *s;

  double T, T_reduced;

  T = 0;
  T_reduced = 0;

  if (solid == -1) {
    for (int isolid = 0; isolid < domain->solids.size(); isolid++) {
      s = domain->solids[isolid];

      for (int ip = 0; ip < s->np_local; ip++) {
        if (s->mask[ip] & groupbit) {
          T += s->T[ip];
        }
      }
    }
  } else {
    s = domain->solids[solid];

    for (int ip = 0; ip < s->np_local; ip++) {
      if (s->mask[ip] & groupbit) {
        T += s->T[ip];
      }
    }
  }

  // Reduce T:
  MPI_Allreduce(&T, &T_reduced, 1, MPI_DOUBLE, MPI_SUM, universe->uworld);
  (*input->vars)[id] = Var(id, T_reduced/group->n_tot(igroup));
}
