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

#include <compute_average_stress.h>
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

ComputeAverageStress::ComputeAverageStress(MPM *mpm, vector<string> args)
    : Compute(mpm, args) {

  if (args.size() < 3) {
    error->all(FLERR, "Error: too few arguments for compute_average_stress: "
                      "requires at least 3 arguments. " +
                          to_string(args.size()) + " received.\n");
  }

  if (group->pon[igroup].compare("particles") != 0 &&
      group->pon[igroup].compare("all") != 0) {
    error->all(FLERR,
               "compute_average_stress needs to be given a group of particles" +
                   group->pon[igroup] + ", " + args[2] + " is a group of " +
                   group->pon[igroup] + ".\n");
  }

  if (universe->me == 0)
    cout << "Creating new compute ComputeAverageStress with ID: " << args[0] << endl;
  id = args[0];

  (*input->vars)[id + "_11"]=Var(id + "_11", 0);
  (*input->vars)[id + "_22"]=Var(id + "_22", 0);
  (*input->vars)[id + "_33"]=Var(id + "_33", 0);
  (*input->vars)[id + "_12"]=Var(id + "_12", 0);
  (*input->vars)[id + "_13"]=Var(id + "_13", 0);
  (*input->vars)[id + "_23"]=Var(id + "_23", 0);
}

ComputeAverageStress::~ComputeAverageStress() {}

void ComputeAverageStress::init() {}

void ComputeAverageStress::setup() {}

void ComputeAverageStress::compute_value() {
  if (update->ntimestep != output->next && update->ntimestep != update->nsteps)
    return;
  // cout << "In ComputeAverageStress::post_particles_to_grid()\n";

  // Go through all the nodes in the group and set b to the right value:

  int solid = group->solid[igroup];

  Solid *s;

  double sigma[6] = {0, 0, 0, 0, 0, 0};
  double sigma_reduced[6] = {0, 0, 0, 0, 0, 0};

  if (solid == -1) {
    for (int isolid = 0; isolid < domain->solids.size(); isolid++) {
      s = domain->solids[isolid];

      for (int ip = 0; ip < s->np_local; ip++) {
        if (s->mask[ip] & groupbit) {
          sigma[0] += s->sigma[ip](0, 0);
          sigma[1] += s->sigma[ip](1, 1);
          sigma[2] += s->sigma[ip](2, 2);
          sigma[3] += s->sigma[ip](1, 2);
          sigma[4] += s->sigma[ip](0, 2);
          sigma[5] += s->sigma[ip](0, 1);
        }
      }
    }
  } else {
    s = domain->solids[solid];

    for (int ip = 0; ip < s->np_local; ip++) {
      if (s->mask[ip] & groupbit) {
	sigma[0] += s->sigma[ip](0, 0);
	sigma[1] += s->sigma[ip](1, 1);
	sigma[2] += s->sigma[ip](2, 2);
	sigma[3] += s->sigma[ip](1, 2);
	sigma[4] += s->sigma[ip](0, 2);
	sigma[5] += s->sigma[ip](0, 1);
      }
    }
  }

  // Reduce Stress:
  MPI_Allreduce(sigma, sigma_reduced, 6, MPI_DOUBLE, MPI_SUM, universe->uworld);

  (*input->vars)[id + "_11"]=Var(id + "_11", sigma_reduced[0]/group->n_tot(igroup));
  (*input->vars)[id + "_22"]=Var(id + "_22", sigma_reduced[1]/group->n_tot(igroup));
  (*input->vars)[id + "_33"]=Var(id + "_33", sigma_reduced[2]/group->n_tot(igroup));
  (*input->vars)[id + "_12"]=Var(id + "_12", sigma_reduced[3]/group->n_tot(igroup));
  (*input->vars)[id + "_02"]=Var(id + "_02", sigma_reduced[4]/group->n_tot(igroup));
  (*input->vars)[id + "_01"]=Var(id + "_01", sigma_reduced[5]/group->n_tot(igroup));
}
