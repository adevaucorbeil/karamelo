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
#include <matrix.h>
#include <iostream>
#include <string>
#include <vector>
using namespace std;
using namespace FixConst;


FixTemperatureNodes::FixTemperatureNodes(MPM *mpm, vector<string> args) : Fix(mpm, args)
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

  Tvalue = input->parsev(args[3]);

  string previous = args[3];

  // Replace "time" by "time - dt" in the x argument:
  previous = SpecialFunc::replace_all(input->parsev(previous).str(), "time", "(time - dt)");
  Tprevvalue = input->parsev(previous);
}

FixTemperatureNodes::~FixTemperatureNodes()
{
}

void FixTemperatureNodes::init()
{
}

void FixTemperatureNodes::setup()
{
}

void FixTemperatureNodes::setmask() {
  mask = 0;
  mask |= POST_UPDATE_GRID_STATE;
  mask |= POST_VELOCITIES_TO_GRID;
}


void FixTemperatureNodes::post_update_grid_state() {
  // cout << "In FixTemperatureNodes::post_update_grid_state()" << endl;

  // Go through all the nodes in the group and set v_update to the right value:
  double T;
  double T_old;

  T = Tvalue.result(mpm);
  T_old = Tprevvalue.result(mpm);

  int solid = group->solid[igroup];
  Grid *g;

  double DT;

  if (solid == -1) {
    for (int isolid = 0; isolid < domain->solids.size(); isolid++) {
      g = domain->solids[isolid]->grid;

      for (int ip = 0; ip < g->nnodes_local + g->nnodes_ghost; ip++) {
	if (g->mask[ip] & groupbit) {
	  DT = 0;
	  DT = T - g->T_update[ip];
	  g->T_update[ip] = T;
	  g->T[ip] = T_old;
	}
      }
    }
  } else {

    g = domain->solids[solid]->grid;

    for (int ip = 0; ip < g->nnodes_local + g->nnodes_ghost; ip++) {
      if (g->mask[ip] & groupbit) {
	DT = 0;
	DT = T - g->T_update[ip];
	g->T_update[ip] = T;
	g->T[ip] = T_old;
      }
    }
  }

  // // Reduce ftot:
  // MPI_Allreduce(ftot.data(), ftot_reduced.data(), 3, MPI_DOUBLE, MPI_SUM,
  //               universe->uworld);

  // (*input->vars)[id + "_x"] = Var(id + "_x", ftot_reduced[0]);
  // (*input->vars)[id + "_y"] = Var(id + "_y", ftot_reduced[1]);
  // (*input->vars)[id + "_z"] = Var(id + "_z", ftot_reduced[2]);
}

void FixTemperatureNodes::post_velocities_to_grid() {
  // cout << "In FixTemperatureNodes::post_velocities_to_grid()" << endl;

  // Go through all the particles in the group and set v to the right value:
  double T = input->parsev(args[3]).result(mpm);

  int solid = group->solid[igroup];
  Grid *g;

  if (solid == -1) {
    for (int isolid = 0; isolid < domain->solids.size(); isolid++) {
      g = domain->solids[isolid]->grid;

      for (int ip = 0; ip < g->nnodes_local + g->nnodes_ghost; ip++) {
	if (g->mask[ip] & groupbit) {
	  g->T[ip] = T;
	}
      }
    }
  } else {
    g = domain->solids[solid]->grid;

    for (int ip = 0; ip < g->nnodes_local + g->nnodes_ghost; ip++) {
      if (g->mask[ip] & groupbit) {
	g->T[ip] = T;
      }
    }
  }
}

void FixTemperatureNodes::write_restart(ofstream *of) {
  Tvalue.write_to_restart(of);
  Tprevvalue.write_to_restart(of);
}

void FixTemperatureNodes::read_restart(ifstream *ifr) {
  Tvalue.read_from_restart(ifr);
  Tprevvalue.read_from_restart(ifr);
}
