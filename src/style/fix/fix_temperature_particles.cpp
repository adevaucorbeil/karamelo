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

#include <fix_temperature_particles.h>
#include <domain.h>
#include <error.h>
#include <group.h>
#include <input.h>
#include <special_functions.h>
#include <universe.h>
#include <update.h>


using namespace std;
using namespace FixConst;



FixTemperatureParticles::FixTemperatureParticles(MPM *mpm, vector<string> args):
  Fix(mpm, args, INITIAL_INTEGRATE | POST_ADVANCE_PARTICLES)
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
    error->all(FLERR, "Error: too few arguments for fix_temperature_nodes.\n" +
                          usage);
  }

  if (group->pon[igroup].compare("particles") !=0 ) {
    error->one(FLERR, "fix_temperature_nodes needs to be given a group of nodes" +
                          group->pon[igroup] + ", " + args[2] +
                          " is a group of " + group->pon[igroup] + ".\n");
  }
  if (universe->me == 0) {
    cout << "Creating new fix FixTemperatureParticles with ID: " << args[0] << endl;
  }
  id = args[0];


  string time = "time";


  Tvalue = input->parsev(args[3]);

  string previous = args[3];

  // Replace "time" by "time - dt" in the x argument:
  previous = SpecialFunc::replace_all(input->parsev(previous).str(), "time", "(time - dt)");
  Tprevvalue = input->parsev(previous);
}

void FixTemperatureParticles::prepare()
{
  Tprevvalue.result(mpm);
  Tvalue    .result(mpm);
}

void FixTemperatureParticles::reduce()
{
  // Reduce ftot:
  //MPI_Allreduce(ftot.elements, ftot_reduced.elements, 3, MPI_DOUBLE, MPI_SUM,
  //              universe->uworld);

  // (*input->vars)[id + "_x"] = Var(id + "_x", ftot_reduced[0]);
  // (*input->vars)[id + "_y"] = Var(id + "_y", ftot_reduced[1]);
  // (*input->vars)[id + "_z"] = Var(id + "_z", ftot_reduced[2]);
}

void FixTemperatureParticles::initial_integrate(Solid &solid, int ip)
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

  //solid.T_update.at(ip) = Tvalue.result(mpm);
  solid.T.at(ip) = Tprevvalue.result(mpm, true);
  // cout << "v_update for " << n << " particles from solid " << domain->solids[solid]->id << " set." << endl;
}

void FixTemperatureParticles::post_advance_particles(Solid &solid, int ip)
{
  // Go through all the particles in the group and set v to the right value:
  if (!(solid.mask.at(ip) & groupbit))
    return;

  (*input->vars)["x" ] = Var("x",  solid.x .at(ip)[0]);
  (*input->vars)["y" ] = Var("y",  solid.x .at(ip)[1]);
  (*input->vars)["z" ] = Var("z",  solid.x .at(ip)[2]);
  (*input->vars)["x0"] = Var("x0", solid.x0.at(ip)[0]);
  (*input->vars)["y0"] = Var("y0", solid.x0.at(ip)[1]);
  (*input->vars)["z0"] = Var("z0", solid.x0.at(ip)[2]);

  solid.T.at(ip) = Tvalue.result(mpm, true);
  // cout << "v for " << n << " particles from solid " <<
  // domain->solids[solid]->id << " set." << endl;
}

void FixTemperatureParticles::write_restart(ofstream *of) {
  Tvalue.write_to_restart(of);
  Tprevvalue.write_to_restart(of);
}

void FixTemperatureParticles::read_restart(ifstream *ifr) {
  Tvalue.read_from_restart(ifr);
  Tprevvalue.read_from_restart(ifr);
}
