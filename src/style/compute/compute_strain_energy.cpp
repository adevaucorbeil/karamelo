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

#include <iostream>
#include <vector>
#include <string>
#include <matrix.h>
#include <compute_strain_energy.h>
#include <input.h>
#include <group.h>
#include <domain.h>
#include <input.h>
#include <update.h>
#include <output.h>
#include <math_special.h>
#include <universe.h>
#include <error.h>
#include <solid.h>

using namespace std;
using namespace MathSpecial;


ComputeStrainEnergy::ComputeStrainEnergy(MPM *mpm, vector<string> args) : Compute(mpm, args)
{
  if (args.size() < 3) {
    error->all(FLERR,"Error: too few arguments for compute_strain_energy: requires at least 3 arguments. " + to_string(args.size()) + " received.\n");
  }

  if (group->pon[igroup].compare("particles") !=0 && group->pon[igroup].compare("all") !=0) {
    error->all(FLERR, "compute_strain_energy needs to be given a group of particles" + group->pon[igroup] + ", " + args[2] + " is a group of " + group->pon[igroup] + ".\n");
  }

  
  if (universe->me == 0) 
    cout << "Creating new compute ComputeStrainEnergy with ID: " << args[0] << endl;
  id = args[0];

  (*input->vars)[id]=Var(id, 0);
}

void ComputeStrainEnergy::compute_value(Solid &solid) {

  double Es, Es_reduced;

  Es = 0;
  Es_reduced = 0;

  Kokkos::View<Matrix3d*> sigma = solid.sigma;
  Kokkos::View<Matrix3d*> strain_el = solid.strain_el;
  Kokkos::View<double*> vol = solid.vol;
  Kokkos::View<int*> mask = solid.mask;

  int groupbit = this->groupbit;

  if (update->ntimestep == output->next ||
      update->ntimestep == update->nsteps)

    Kokkos::parallel_reduce("ComputeStrainEnergy::compute_value", solid.np_local,
			    KOKKOS_LAMBDA(const int &ip, double &lEs) {
			      if (mask[ip] & groupbit)
				lEs += 0.5*vol[ip]*(sigma[ip](0,0)*strain_el[ip](0,0)
						    + sigma[ip](0,1)*strain_el[ip](0,1)
						    + sigma[ip](0,2)*strain_el[ip](0,2)
						    + sigma[ip](1,0)*strain_el[ip](1,0)
						    + sigma[ip](1,1)*strain_el[ip](1,1)
						    + sigma[ip](1,2)*strain_el[ip](1,2)
						    + sigma[ip](2,0)*strain_el[ip](2,0)
						    + sigma[ip](2,1)*strain_el[ip](2,1)
						    + sigma[ip](2,2)*strain_el[ip](2,2));
			    },Es);


  // Reduce Es:
  MPI_Allreduce(&Es,&Es_reduced,1,MPI_DOUBLE,MPI_SUM,universe->uworld);
  (*input->vars)[id]=Var(id, Es_reduced);
}
