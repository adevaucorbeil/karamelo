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

#include <fix_temperature_nodes.h>
#include <domain.h>
#include <error.h>
#include <grid.h>
#include <group.h>
#include <input.h>
#include <special_functions.h>
#include <universe.h>
#include <update.h>


using namespace std;
using namespace FixConst;


FixTemperatureNodes::FixTemperatureNodes(MPM *mpm, vector<string> args):
  Fix(mpm, args, POST_UPDATE_GRID_STATE | POST_VELOCITIES_TO_GRID)
{
  if (args.size() < Nargs) {
    error->all(FLERR, "Error: not enough arguments.\n" + usage);
  }

  if (group->pon[igroup].compare("nodes") !=0 ) {
    error->one(FLERR, "fix_temperature_nodes needs to be given a group of nodes" + group->pon[igroup] + ", " + args[2] + " is a group of "+ group->pon[igroup] + ".\n");
  }
  if (universe->me == 0) {
    cout << "Creating new fix FixTemperatureNodes with ID: " << args[0] << endl;
  }
  id = args[0];

  string time = "time";

  Targ = args[3];
  Tvalue = input->parsev(args[3]);

  string previous = args[3];

  // Replace "time" by "time - dt" in the x argument:
  previous = SpecialFunc::replace_all(input->parsev(previous).str(), "time", "(time - dt)");
  Tprevvalue = input->parsev(previous);
}

void FixTemperatureNodes::prepare()
{
  Tvalue    .result(mpm);
  Tprevvalue.result(mpm);
}

void FixTemperatureNodes::reduce()
{
  // // Reduce ftot:
  // MPI_Allreduce(ftot.elements, ftot_reduced.elements, 3, MPI_DOUBLE, MPI_SUM,
  //               universe->uworld);

  // (*input->vars)[id + "_x"] = Var(id + "_x", ftot_reduced[0]);
  // (*input->vars)[id + "_y"] = Var(id + "_y", ftot_reduced[1]);
  // (*input->vars)[id + "_z"] = Var(id + "_z", ftot_reduced[2]);
}

void FixTemperatureNodes::post_update_grid_state(Grid &grid, int in)
{
  // cout << "In FixTemperatureNodes::post_update_grid_state()" << endl;

  // Go through all the nodes in the group and set v_update to the right value:
  if (!(grid.mask[in] & groupbit))
    return;

  grid.T_update[in] = Tvalue    .result(mpm, true);
  grid.T       [in] = Tprevvalue.result(mpm, true);
}

void FixTemperatureNodes::post_velocities_to_grid(Grid &grid, int in)
{
  // cout << "In FixTemperatureNodes::post_velocities_to_grid()" << endl;

  // Go through all the particles in the group and set v to the right value:
  if (grid.mask[in] & groupbit)
    grid.T[in] = input->parsev(Targ).result(mpm);
}

void FixTemperatureNodes::write_restart(ofstream *of) {
  Tvalue.write_to_restart(of);
  Tprevvalue.write_to_restart(of);
}

void FixTemperatureNodes::read_restart(ifstream *ifr) {
  Tvalue.read_from_restart(ifr);
  Tprevvalue.read_from_restart(ifr);
}
