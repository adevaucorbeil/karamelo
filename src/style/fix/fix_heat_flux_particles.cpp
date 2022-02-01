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

  if (args[3] != "NULL") {
    q = input->parsev(args[3]);
  }
}

void FixHeatFluxParticles::prepare()
{
  q.result(mpm);

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

void FixHeatFluxParticles::initial_integrate(Solid &solid, int ip)
{
  // Go through all the particles in the group and set v_update to the right value:
  if (!(solid.mask.at(ip) & groupbit))
    return;

  (*input->vars)["x" ] = Var("x",  solid.x .at(ip)[0]);
  (*input->vars)["y" ] = Var("y",  solid.x .at(ip)[1]);
  (*input->vars)["z" ] = Var("z",  solid.x .at(ip)[2]);
  (*input->vars)["x0"] = Var("x0", solid.x0.at(ip)[0]);
  (*input->vars)["y0"] = Var("y0", solid.x0.at(ip)[1]);
  (*input->vars)["z0"] = Var("z0", solid.x0.at(ip)[2]);

  double Ap;
  if (domain->dimension == 1)
    Ap = 1;
  else if (domain->dimension == 2)
    Ap = sqrt(solid.vol.at(ip));
  else
    Ap = pow(solid.vol.at(ip), 2/3);

  double qtemp = q.result(mpm, true);
  solid.gamma.at(ip) += Ap*qtemp*solid.mat->invcp;
  qtot += qtemp;
}

void FixHeatFluxParticles::write_restart(ofstream *of) {
  q.write_to_restart(of);
}

void FixHeatFluxParticles::read_restart(ifstream *ifr) {
  q.read_from_restart(ifr);
}
