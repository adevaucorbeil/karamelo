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

#include <fix_initial_velocity_nodes.h>
#include <domain.h>
#include <error.h>
#include <grid.h>
#include <group.h>
#include <input.h>
#include <universe.h>
#include <update.h>


using namespace std;
using namespace FixConst;


FixInitialVelocityNodes::FixInitialVelocityNodes(MPM *mpm, vector<string> args):
  Fix(mpm, args, POST_UPDATE_GRID_STATE | POST_VELOCITIES_TO_GRID)
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

  if (domain->dimension == 3 && args.size()<6) {
    error->all(FLERR,"Error: too few arguments for fix_initial_velocity_nodes: requires at least 6 arguments. " + to_string(args.size()) + " received.\n");
  } else if (domain->dimension == 2 && args.size()<5) {
    error->all(FLERR,"Error: too few arguments for fix_initial_velocity_nodes: requires at least 5 arguments. " + to_string(args.size()) + " received.\n");
  } else if (domain->dimension == 1 && args.size()<4) {
    error->all(FLERR,"Error: too few arguments for fix_initial_velocity_nodes: requires at least 4 arguments. " + to_string(args.size()) + " received.\n");
  }

  if (group->pon[igroup].compare("nodes") !=0 && group->pon[igroup].compare("all") !=0) {
    error->all(FLERR,"fix_initial_velocity_nodes needs to be given a group of nodes" + group->pon[igroup] + ", " + args[2] + " is a group of " + group->pon[igroup] + ".\n");
  }
  if (universe->me == 0) {
    cout << "Creating new fix FixInitialVelocityNodes with ID: " << args[0] << endl;
  }
  id = args[0];

  xset = yset = zset = false;

  if (args[3] != "NULL") {
    xvalue = input->parsev(args[3]);
    xset = true;
  }

  if (domain->dimension >= 2) {
    if (args[4] != "NULL") {
      yvalue = input->parsev(args[4]);
      yset = true;
    }
  }

  if (domain->dimension == 3) {
    if (args[5] != "NULL") {
      zvalue = input->parsev(args[5]);
      zset = true;
    }
  }
}

void FixInitialVelocityNodes::prepare()
{
  xvalue.result(mpm);
  yvalue.result(mpm);
  zvalue.result(mpm);
}

void FixInitialVelocityNodes::post_update_grid_state(Grid &grid, int in)
{
  // cout << "In FixInitialVelocityNodes::post_update_grid_state()" << endl;

  // Go through all the nodes in the group and set v_update to the right value:
  if (update->ntimestep != 1 || !(grid.mask.at(in) & groupbit))
    return;
  
  (*input->vars)["x0"] = Var("x0", grid.x0.at(in)[0]);
  (*input->vars)["y0"] = Var("y0", grid.x0.at(in)[1]);
  (*input->vars)["z0"] = Var("z0", grid.x0.at(in)[2]);

  if (xset)
	grid.v_update.at(in)[0] = xvalue.result(mpm, true);
  if (yset)
	grid.v_update.at(in)[1] = yvalue.result(mpm, true);
  if (zset)
	grid.v_update.at(in)[2] = zvalue.result(mpm, true);
  // cout << "grid.v_update for " << n << " nodes from solid " << domain->solids[solid]->id << " set." << endl;
}

void FixInitialVelocityNodes::post_velocities_to_grid(Grid &grid, int in)
{
  // cout << "In FixInitialVelocityNodes::post_velocities_to_grid()" << endl;

  // Go through all the particles in the group and set v to the right value:
  if (update->ntimestep != 1 || !(grid.mask.at(in) & groupbit))
    return;
  
  (*input->vars)["x0"] = Var("x0", grid.x0.at(in)[0]);
  (*input->vars)["y0"] = Var("y0", grid.x0.at(in)[1]);
  (*input->vars)["z0"] = Var("z0", grid.x0.at(in)[2]);

  if (xset)
	grid.v.at(in)[0] = xvalue.result(mpm, true);
  if (yset)
	grid.v.at(in)[1] = yvalue.result(mpm, true);
  if (zset)
	grid.v.at(in)[2] = zvalue.result(mpm, true);
  // cout << "v for " << n << " nodes from solid " << domain->solids[solid]->id << " set." << endl;
}
