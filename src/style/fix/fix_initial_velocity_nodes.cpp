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
#include <expression_operation.h>
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

  if (args.size() < Nargs.find(domain->dimension)->second) {
    error->all(FLERR, "Error: too few arguments for fix_initial_velocity_nodes.\n" +
                          usage.find(domain->dimension)->second);
  }

  if (group->pon[igroup].compare("nodes") !=0 && group->pon[igroup].compare("all") !=0) {
    error->all(FLERR,"fix_initial_velocity_nodes needs to be given a group of nodes" + group->pon[igroup] + ", " + args[2] + " is a group of " + group->pon[igroup] + ".\n");
  }
  if (universe->me == 0) {
    cout << "Creating new fix FixInitialVelocityNodes with ID: " << args[0] << endl;
  }
  id = args[0];

  v[0] = nullptr;
  v[1] = nullptr;
  v[2] = nullptr;

  
  for (int dim = 0; dim < domain->dimension; dim++) {
    if (args[3 + dim] != "NULL") {
      input->parsev(args[3 + dim]);
      v[dim] = &input->expressions[args[3 + dim]];
    }
  }
}

void FixInitialVelocityNodes::prepare()
{
}

void FixInitialVelocityNodes::post_update_grid_state(Grid &grid)
{
  if (update->ntimestep != 1)
    return;

  // Go through all the nodes in the group and set v_update to the right value:
  for (int i = 0; i < 3; i++) {
    if (v[i])
      v[i]->evaluate(grid);
  }

  int groupbit = this->groupbit;
  Kokkos::View<int*> mask = grid.mask;
  Kokkos::View<Vector3d*> gv_update = grid.v_update;

  for (int i = 0; i < 3; i++)
    if (v[i])
    {
      Kokkos::View<float **> v_i = v[i]->registers;

      Kokkos::parallel_for("FixInitialVelocityNodes::post_update_grid_state", grid.nnodes_local + grid.nnodes_ghost,
      KOKKOS_LAMBDA(const int &in)
      {
        if (!(mask[in] & groupbit))
          return;
	gv_update[in][i] = v_i(0, in);
      });
    }
}

void FixInitialVelocityNodes::post_velocities_to_grid(Grid &grid)
{
  if (update->ntimestep != 1)
    return;
  // cout << "In FixInitialVelocityNodes::post_velocities_to_grid()" << endl;

  // Go through all the particles in the group and set v to the right value:
  for (int i = 0; i < 3; i++) {
    if (v[i])
      v[i]->evaluate(grid);
  }

  int groupbit = this->groupbit;
  Kokkos::View<int*> mask = grid.mask;
  Kokkos::View<Vector3d*> gv = grid.v;

  for (int i = 0; i < 3; i++)
    if (v[i])
    {
      Kokkos::View<float **> v_i = v[i]->registers;

      Kokkos::parallel_for("FixInitialVelocityNodes::post_update_grid_state", grid.nnodes_local + grid.nnodes_ghost,
      KOKKOS_LAMBDA(const int &in)
      {
        if (!(mask[in] & groupbit))
          return;
	gv[in][i] = v_i(0, in);
      });
    }
}
