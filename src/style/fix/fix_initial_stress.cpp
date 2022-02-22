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

#include <fix_initial_stress.h>
#include <domain.h>
#include <error.h>
#include <group.h>
#include <input.h>
#include <solid.h>
#include <universe.h>
#include <update.h>


using namespace std;
using namespace FixConst;


FixInitialStress::FixInitialStress(MPM *mpm, vector<string> args):
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
    error->all(FLERR, "fix_initial_stress needs to be given a group of particles" + group->pon[igroup] + ", " + args[2] + " is a group of " + group->pon[igroup] + ".\n");
  }
  if (universe->me == 0) {
    cout << "Creating new fix FixInitialStress with ID: " << args[0] << endl;
  }
  id = args[0];

  for (int i = 0; i < 6; i++) {
    if (args[i + 3] != "NULL") {
      s_value[i] = input->parsev(args[i + 3]);
      s_set[i] = true;
    } else {
      s_set[i] = false;
    }
  }
}

void FixInitialStress::prepare()
{
  for (Var &s_value: s_value)
    s_value.result(mpm);
}

void FixInitialStress::initial_integrate(Solid &solid, int ip)
{
  // cout << "In FixInitialStress::initial_integrate()" << endl;

  // Go through all the particles in the group and set v to the right value:
  if (update->ntimestep != 1 || !(solid.mask.at(ip) & groupbit))
    return;

  (*input->vars)["x0"] = Var("x0", solid.x0.at(ip)[0]);
  (*input->vars)["y0"] = Var("y0", solid.x0.at(ip)[1]);
  (*input->vars)["z0"] = Var("z0", solid.x0.at(ip)[2]);
  (*input->vars)["x" ] = Var("x",  solid.x[ip][0]);
  (*input->vars)["y" ] = Var("y",  solid.x[ip][1]);
  (*input->vars)["z" ] = Var("z",  solid.x[ip][2]);
      
  if (s_set[0]) solid.sigma.at(ip)(0,0) = s_value[0].result(mpm, true);
  if (s_set[1]) solid.sigma.at(ip)(1,1) = s_value[1].result(mpm, true);
  if (s_set[2]) solid.sigma.at(ip)(2,2) = s_value[2].result(mpm, true);
  if (s_set[3]) solid.sigma.at(ip)(1,2) = solid.sigma.at(ip)(2,1) = s_value[3].result(mpm, true);
  if (s_set[4]) solid.sigma.at(ip)(0,2) = solid.sigma.at(ip)(2,0) = s_value[4].result(mpm, true);
  if (s_set[5]) solid.sigma.at(ip)(0,1) = solid.sigma.at(ip)(1,0) = s_value[5].result(mpm, true);

  if (update->method_type == "tlmpm" || update->method_type == "tlcpdi")
    solid.vol0PK1.at(ip) = solid.vol0.at(ip)*solid.sigma.at(ip);
}
