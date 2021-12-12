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
#include <Eigen/Eigen>
#include <iostream>
#include <vector>

using namespace std;
using namespace FixConst;
using namespace Eigen;


FixTemperatureParticles::FixTemperatureParticles(MPM *mpm, vector<string> args) : Fix(mpm, args)
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

FixTemperatureParticles::~FixTemperatureParticles()
{
}

void FixTemperatureParticles::init()
{
}

void FixTemperatureParticles::setup()
{
}

void FixTemperatureParticles::setmask() {
  mask = 0;
  mask |= INITIAL_INTEGRATE;
  mask |= POST_ADVANCE_PARTICLES;
}


void FixTemperatureParticles::initial_integrate() {
  // Go through all the particles in the group and set v_update to the right value:
  int solid = group->solid[igroup];
  Solid *s;

  int n = 0;

  if (solid == -1) {
    for (int isolid = 0; isolid < domain->solids.size(); isolid++) {
      s = domain->solids[isolid];
      n = 0;

      for (int ip = 0; ip < s->np_local; ip++) {
	if (s->mask[ip] & groupbit) {
	  (*input->vars)["x"] = Var("x", s->x[ip][0]);
	  (*input->vars)["x0"] = Var("x0", s->x0[ip][0]);
	  (*input->vars)["y"] = Var("y", s->x[ip][1]);
	  (*input->vars)["y0"] = Var("y0", s->x0[ip][1]);
	  (*input->vars)["z"] = Var("z", s->x[ip][2]);
	  (*input->vars)["z0"] = Var("z0", s->x0[ip][2]);

	  //s->T_update[ip] = Tvalue.result(mpm);
	  s->T[ip] = Tprevvalue.result(mpm);
	  n++;
	}
      }
      // cout << "v_update for " << n << " particles from solid " << domain->solids[isolid]->id << " set." << endl;
    }
  } else {
    s = domain->solids[solid];

    for (int ip = 0; ip < s->np_local; ip++) {
      if (s->mask[ip] & groupbit) {
	(*input->vars)["x"] = Var("x", s->x[ip][0]);
	(*input->vars)["x0"] = Var("x0", s->x0[ip][0]);
	(*input->vars)["y"] = Var("y", s->x[ip][1]);
	(*input->vars)["y0"] = Var("y0", s->x0[ip][1]);
	(*input->vars)["z"] = Var("z", s->x[ip][2]);
	(*input->vars)["z0"] = Var("z0", s->x0[ip][2]);


	//s->T_update[ip] = Tvalue.result(mpm);
	s->T[ip] = Tprevvalue.result(mpm);
	n++;
      }
    }
    // cout << "v_update for " << n << " particles from solid " << domain->solids[solid]->id << " set." << endl;
  }
}

void FixTemperatureParticles::post_advance_particles() {
  // Go through all the particles in the group and set v to the right value:
  int solid = group->solid[igroup];
  Solid *s;
  // Eigen::Vector3d ftot, ftot_reduced;

  int n = 0;
  //ftot.setZero();
  double inv_dt = 1.0/update->dt;

  if (solid == -1) {
    for (int isolid = 0; isolid < domain->solids.size(); isolid++) {
      s = domain->solids[isolid];
      n = 0;

      for (int ip = 0; ip < s->np_local; ip++) {
        if (s->mask[ip] & groupbit) {
	  (*input->vars)["x"] = Var("x", s->x[ip][0]);
	  (*input->vars)["x0"] = Var("x0", s->x0[ip][0]);
	  (*input->vars)["y"] = Var("y", s->x[ip][1]);
	  (*input->vars)["y0"] = Var("y0", s->x0[ip][1]);
	  (*input->vars)["z"] = Var("z", s->x[ip][2]);
	  (*input->vars)["z0"] = Var("z0", s->x0[ip][2]);

	  s->T[ip] = Tvalue.result(mpm);
          n++;
        }
      }
      // cout << "v for " << n << " particles from solid " <<
      // domain->solids[isolid]->id << " set." << endl;
    }
  } else {
    s = domain->solids[solid];
    n = 0;
    for (int ip = 0; ip < s->np_local; ip++) {
      if (s->mask[ip] & groupbit) {
	(*input->vars)["x"] = Var("x", s->x[ip][0]);
	(*input->vars)["x0"] = Var("x0", s->x0[ip][0]);
	(*input->vars)["y"] = Var("y", s->x[ip][1]);
	(*input->vars)["y0"] = Var("y0", s->x0[ip][1]);
	(*input->vars)["z"] = Var("z", s->x[ip][2]);
	(*input->vars)["z0"] = Var("z0", s->x0[ip][2]);

	s->T[ip] = Tvalue.result(mpm);
        n++;
      }
    }
    // cout << "v for " << n << " particles from solid " <<
    // domain->solids[solid]->id << " set." << endl;
  }

  // Reduce ftot:
  //MPI_Allreduce(ftot.data(), ftot_reduced.data(), 3, MPI_DOUBLE, MPI_SUM,
  //              universe->uworld);

  // (*input->vars)[id + "_x"] = Var(id + "_x", ftot_reduced[0]);
  // (*input->vars)[id + "_y"] = Var(id + "_y", ftot_reduced[1]);
  // (*input->vars)[id + "_z"] = Var(id + "_z", ftot_reduced[2]);
}

void FixTemperatureParticles::write_restart(ofstream *of) {
  Tvalue.write_to_restart(of);
  Tprevvalue.write_to_restart(of);
}

void FixTemperatureParticles::read_restart(ifstream *ifr) {
  Tvalue.read_from_restart(ifr);
  Tprevvalue.read_from_restart(ifr);
}
