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
#include <expression_operation.h>


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
  string previous = args[3];

  input->parsev(args[3]);
  Tvalue = &input->expressions[args[3]];


  // Replace "time" by "time - dt" in the x argument:
  previous = SpecialFunc::replace_all(input->parsev(previous).str(), "time", "(time - dt)");
  input->parsev(previous);
  Tprevvalue = &input->expressions[previous];
}

void FixTemperatureNodes::post_update_grid_state(Grid &grid)
{
  // Update the temperatures:
  Tvalue->evaluate(grid);
  Tprevvalue->evaluate(grid);

  int groupbit = this->groupbit;
  Kokkos::View<int*> mask = grid.mask;
  Kokkos::View<float*> T = grid.T, T_update = grid.T_update;

  Kokkos::View<float **> Tv = Tvalue->registers;
  Kokkos::View<float **> Tpv = Tprevvalue->registers;

  Kokkos::parallel_for("FixTemperatureNodes::post_update_grid_state", grid.nnodes_local + grid.nnodes_ghost,
		       KOKKOS_LAMBDA(const int &in)
      {
        if (!(mask[in] & groupbit))
          return;

        T_update[in] = Tv(0, in);
        T[in]        = Tpv(0, in);
      });
}

void FixTemperatureNodes::post_velocities_to_grid(Grid &grid)
{
  // Update the temperatures:
  Tvalue->evaluate(grid);

  int groupbit = this->groupbit;
  Kokkos::View<int*> mask = grid.mask;
  Kokkos::View<float*> T = grid.T;

  Kokkos::View<float **> Tv = Tvalue->registers;

  Kokkos::parallel_for("FixTemperatureNodes::post_velocities_to_grid", grid.nnodes_local + grid.nnodes_ghost,
		       KOKKOS_LAMBDA(const int &in)
      {
        if (!(mask[in] & groupbit))
          return;

        T[in] = Tv(0, in);
      });
}

void FixTemperatureNodes::write_restart(ofstream *of) {
  // Tvalue.write_to_restart(of);
  // Tprevvalue.write_to_restart(of);
}

void FixTemperatureNodes::read_restart(ifstream *ifr) {
  // Tvalue.read_from_restart(ifr);
  // Tprevvalue.read_from_restart(ifr);
}
