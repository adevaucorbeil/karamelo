/* ----------------------------------------------------------------------
 *
 *                    ***       Karamelo       ***
 *               Parallel Material Point Method Simulator
 * 
 * Copyright (2020) Alban de Vaucorbeil, alban.devaucorbeil@deakin.edu.au
 * Institute for Frontier Materials, Deakin University
 * Geelong VIC 3216, Australia

 * This software is distributed under the GNU General Public License.
 *
 * ----------------------------------------------------------------------- */

#include <fix_heat_flux_particles.h>
#include <domain.h>
#include <error.h>
#include <group.h>
#include <input.h>
#include <special_functions.h>
#include <universe.h>
#include <update.h>
#include <expression_operation.h>

using namespace std;
using namespace FixConst;


FixHeatFluxParticles::FixHeatFluxParticles(MPM *mpm, vector<string> args):
  Fix(mpm, args, INITIAL_INTEGRATE)
{
  if (args.size() < 3) {
    error->all(FLERR, "Error: not enough arguments.\n");
  }

  if (args[2].compare("restart") ==
      0) { // If the keyword restart, we are expecting to have read_restart()
           // launched right after.
    igroup = stoi(args[3]);
    if (igroup == -1 && universe->me == 0) {
      cout << "Could not find group number " << args[3] << endl;
    }
    groupbit = group->bitmask[igroup];
    return;
  }

  if (args.size() < Nargs) {
    error->all(FLERR, "Error: too few arguments for fix_heat_flux_particles.\n" +
                          usage);
  }

  if (args.size() > Nargs) {
    error->all(FLERR, "Error: too many arguments for fix_heat_flux_particles.\n" +
           usage);
  }

  if (igroup == -1) {
    error->all(FLERR, "Could not find group ID " + args[2] + "\n");
  }

  if (group->pon[igroup] != "particles") {
    error->one(FLERR, "fix_heat_flux_nodes needs to be given a group of particles" +
                          group->pon[igroup] + ", " + args[2] +
                          " is a group of " + group->pon[igroup] + ".\n");
  }
  if (universe->me == 0) {
    cout << "Creating new fix FixHeatFluxParticles with ID: " << args[0] << endl;
  }

  id = args[0];
  
  input->parsev(args[3]);
  q = &input->expressions[args[3]];
}

void FixHeatFluxParticles::prepare()
{
  q->evaluate();
  qtot = 0;
}

void FixHeatFluxParticles::reduce()
{
  double qtot_reduced;

  // Reduce qtot:
  MPI_Allreduce(&qtot, &qtot_reduced, 1, MPI_DOUBLE, MPI_SUM,
                universe->uworld);

  (*input->vars)[id + "_s"] = Var(id + "_s", qtot_reduced);
}

void FixHeatFluxParticles::initial_integrate(Solid &solid)
{
  // Go through all the particles in the group and set v_update to the right value:

  int groupbit = this->groupbit, dimension = domain->dimension;
  Kokkos::View<int*> mask = solid.mask;
  Kokkos::View<double*> T = solid.T, vol = solid.vol, gamma = solid.gamma;
  double invcp = solid.mat->invcp;

  Kokkos::View<double**> q_ = q->registers;

  Kokkos::parallel_reduce("FixVelocityNodes::post_update_grid_state", solid.np_local,
			  KOKKOS_LAMBDA(const int &ip, double &lqtot)
      {
        if (!(mask[ip] & groupbit))
          return;

	double Ap;
	if (dimension == 1)
	  Ap = 1;
	else if (dimension == 2)
	  Ap = Kokkos::Experimental::sqrt(vol[ip]);
	else         
	  Ap = Kokkos::Experimental::pow(vol[ip], 2/3);

	gamma[ip] += Ap * q_(0,ip) * invcp;
	lqtot += q_(0,ip);
      }, qtot);
}

void FixHeatFluxParticles::write_restart(ofstream *of) {
  // q.write_to_restart(of);
}

void FixHeatFluxParticles::read_restart(ifstream *ifr) {
  // q.read_from_restart(ifr);
}
