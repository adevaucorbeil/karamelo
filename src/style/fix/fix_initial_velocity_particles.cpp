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

  if (args.size() < Nargs) {
    error->all(FLERR, "Error: not enough arguments.\n" + usage);
  }

  if (args.size() > Nargs) {
    error->all(FLERR, "Error: too many arguments.\n" + usage);
  }

  if (group->pon[igroup].compare("particles") !=0 && group->pon[igroup].compare("all") !=0) {
    error->all(FLERR, "fix_initial_velocity_particles needs to be given a group of particles" + group->pon[igroup] + ", " + args[2] + " is a group of " + group->pon[igroup] + ".\n");
  }

  id = args[0];

  xset = yset = zset = false;

  if (args[3] != "NULL") {
    xvalue = input->parsev(args[3]);
    xset = true;
  }

  if (args[4] != "NULL") {
    yvalue = input->parsev(args[4]);
    yset = true;
  }

  if (args[5] != "NULL") {
    zvalue = input->parsev(args[5]);
    zset = true;
  }
}

void FixInitialVelocityParticles::prepare()
{
  xvalue.result(mpm);
  yvalue.result(mpm);
  zvalue.result(mpm);
}

void FixInitialVelocityParticles::initial_integrate(Solid &solid, int ip) {
  // cout << "In FixInitialVelocityParticles::initial_integrate()" << endl;

  // Go through all the particles in the group and set v to the right value:
  if (update->ntimestep != 1 || !(solid.mask.at(ip) & groupbit))
    return;

  (*input->vars)["x" ] = Var("x",  solid.x .at(ip)[0]);
  (*input->vars)["y" ] = Var("y",  solid.x .at(ip)[1]);
  (*input->vars)["z" ] = Var("z",  solid.x .at(ip)[2]);
  (*input->vars)["x0"] = Var("x0", solid.x0.at(ip)[0]);
  (*input->vars)["y0"] = Var("y0", solid.x0.at(ip)[1]);
  (*input->vars)["z0"] = Var("z0", solid.x0.at(ip)[2]);

  if (xset)
    solid.v.at(ip)[0] = xvalue.result(mpm, true);
  if (yset)
	solid.v.at(ip)[1] = yvalue.result(mpm, true);
  if (zset)
	solid.v.at(ip)[2] = zvalue.result(mpm, true);
}
