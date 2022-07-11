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

#include <fix_initial_velocity_particles.h>
#include <domain.h>
#include <error.h>
#include <group.h>
#include <input.h>
#include <solid.h>
#include <universe.h>
#include <update.h>
#include <expression_operation.h>


using namespace std;
using namespace FixConst;


FixInitialVelocityParticles::FixInitialVelocityParticles(MPM *mpm, vector<string> args):
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

  if (args.size() < Nargs.find(domain->dimension)->second) {
    error->all(FLERR, "Error: too few arguments for fix_velocity_nodes.\n" +
                          usage.find(domain->dimension)->second);
  }

  if (group->pon[igroup].compare("particles") !=0 && group->pon[igroup].compare("all") !=0) {
    error->all(FLERR, "fix_initial_velocity_particles needs to be given a group of particles" + group->pon[igroup] + ", " + args[2] + " is a group of " + group->pon[igroup] + ".\n");
  }

  if (universe->me == 0) {
    cout << "Creating new fix FixInitialVelocityParticles with ID: " << args[0] << endl;
  }

  id = args[0];

  xset = yset = zset = false;

  if (args[3] != "NULL") {
    input->parsev(args[3]);
    xset = true;
    v[0] = &input->expressions[args[3]];
  }
  else
    v[0] = nullptr;

  if (domain->dimension >= 2 && args[4] != "NULL") {
    input->parsev(args[4]);
    yset = true;
    v[1] = &input->expressions[args[4]];
  }
  else
    v[1] = nullptr;

  if (domain->dimension == 3 && args[5] != "NULL") {
    input->parsev(args[5]);
    zset = true;
    v[2] = &input->expressions[args[5]];
  } else
    v[2] = nullptr;
}

void FixInitialVelocityParticles::initial_integrate(Solid &solid) {
  // cout << "In FixInitialVelocityParticles::initial_integrate()" << endl;

  // Go through all the nodes in the group and set v to the right value:

  if (update->ntimestep != 1)
    return;

  for (int i = 0; i < 3; i++)
    if (v[i])
      v[i]->evaluate(solid);


  int groupbit = this->groupbit;
  Kokkos::View<int*> mask = solid.mask;
  Kokkos::View<Vector3d*> sv = solid.v;

  for (int i = 0; i < 3; i++)
    if (v[i])
      {
	Kokkos::View<float **> v_i = v[i]->registers;

	Kokkos::parallel_for("FixInitialVelocityParticles::initial_integrate", solid.np_local,
			     KOKKOS_LAMBDA(const int &ip)
			     {
			       if (mask[ip] & groupbit)
				 sv[ip][i] = v_i(0, ip);
			     });
      }
}
