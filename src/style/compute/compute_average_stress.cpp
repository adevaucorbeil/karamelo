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
#include <iostream>
#include <math_special.h>
#include <output.h>
#include <solid.h>
#include <string>
#include <universe.h>
#include <update.h>
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

void ComputeAverageStress::compute_value(Solid &solid) {

  double s0, s1, s2, s3, s4, s5;

  Kokkos::View<Matrix3d*> sigma = solid.sigma;
  Kokkos::View<int*> mask = solid.mask;

  int groupbit = this->groupbit;
  int nsteps = update->nsteps;
  bigint next = output->next;
  double ntimestep = update->ntimestep;

  Kokkos::parallel_reduce("ComputeAverageStress::compute_value", solid.np_local,
			  KOKKOS_LAMBDA(const int &ip,
					double &ls0, double &ls1, double &ls2,
					double &ls3, double &ls4, double &ls5) {
			 if ((ntimestep != next &&
			      ntimestep != nsteps) ||
			     !(mask[ip] & groupbit))
			   return;

			 ls0 += sigma[ip](0, 0);
			 ls1 += sigma[ip](1, 1);
			 ls2 += sigma[ip](2, 2);
			 ls3 += sigma[ip](1, 2);
			 ls4 += sigma[ip](0, 2);
			 ls5 += sigma[ip](0, 1);
			  },s0, s1, s2, s3, s4, s5);

  // Reduce Stress:

  double s[6] = {s0, s1, s2, s3, s4, s5};
  double sigma_reduced[6] = {0, 0, 0, 0, 0, 0};
  MPI_Allreduce(s, sigma_reduced, 6, MPI_DOUBLE, MPI_SUM, universe->uworld);

  (*input->vars)[id + "_11"]=Var(id + "_11", sigma_reduced[0]/group->n_tot(igroup));
  (*input->vars)[id + "_22"]=Var(id + "_22", sigma_reduced[1]/group->n_tot(igroup));
  (*input->vars)[id + "_33"]=Var(id + "_33", sigma_reduced[2]/group->n_tot(igroup));
  (*input->vars)[id + "_12"]=Var(id + "_12", sigma_reduced[3]/group->n_tot(igroup));
  (*input->vars)[id + "_02"]=Var(id + "_02", sigma_reduced[4]/group->n_tot(igroup));
  (*input->vars)[id + "_01"]=Var(id + "_01", sigma_reduced[5]/group->n_tot(igroup));
}
