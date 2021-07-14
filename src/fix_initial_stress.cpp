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

#include "fix_initial_stress.h"
#include "domain.h"
#include "error.h"
#include "group.h"
#include "input.h"
#include "solid.h"
#include "universe.h"
#include "update.h"
#include <Eigen/Eigen>
#include <iostream>
#include <string>
#include <vector>

using namespace std;
using namespace FixConst;
using namespace Eigen;

FixInitialStress::FixInitialStress(MPM *mpm, vector<string> args) : Fix(mpm, args)
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
    if (args[i + 3].compare("NULL") != 0) {
      s_value[i] = input->parsev(args[i + 3]);
      s_set[i] = true;
    } else {
      s_set[i] = false;
    }
  }
}

FixInitialStress::~FixInitialStress()
{
}

void FixInitialStress::init()
{
}

void FixInitialStress::setup()
{
}

void FixInitialStress::setmask() {
  mask = 0;
  mask |= INITIAL_INTEGRATE;
}


void FixInitialStress::initial_integrate() {
  if (update->ntimestep !=1) return;
  // cout << "In FixInitialStress::initial_integrate()" << endl;

  // Go through all the particles in the group and set v to the right value:

  int solid = group->solid[igroup];

  Solid *s;

  bool tl;

  if (update->method_type.compare("tlmpm") == 0 ||
      update->method_type.compare("tlcpdi") == 0)
    tl = true;
  else
    tl = false;

  if (solid == -1) {
    for (int isolid = 0; isolid < domain->solids.size(); isolid++) {
      s = domain->solids[isolid];

      for (int ip = 0; ip < s->np_local; ip++) {
	if (s->mask[ip] & groupbit) {
	  (*input->vars)["x0"] = Var("x0", s->x0[ip][0]);
	  (*input->vars)["y0"] = Var("y0", s->x0[ip][1]);
	  (*input->vars)["z0"] = Var("z0", s->x0[ip][2]);
	  (*input->vars)["x"] = Var("x", s->x[ip][0]);
	  (*input->vars)["y"] = Var("y", s->x[ip][1]);
	  (*input->vars)["z"] = Var("z", s->x[ip][2]);
	  
	  if (s_set[0])  s->sigma[ip](0,0) = s_value[0].result(mpm);
	  if (s_set[1])  s->sigma[ip](1,1) = s_value[1].result(mpm);
	  if (s_set[2])  s->sigma[ip](2,2) = s_value[2].result(mpm);
	  if (s_set[3])  s->sigma[ip](1,2) = s->sigma[ip](2,1) = s_value[3].result(mpm);
	  if (s_set[4])  s->sigma[ip](0,2) = s->sigma[ip](2,0) = s_value[4].result(mpm);
	  if (s_set[5])  s->sigma[ip](0,1) = s->sigma[ip](1,0) = s_value[5].result(mpm);

	  if (tl) {
	    s->vol0PK1[ip] = s->vol0[ip] * s->sigma[ip];
	  }
	}
      }
    }
  } else {
    s = domain->solids[solid];

    for (int ip = 0; ip < s->np_local; ip++) {
      if (s->mask[ip] & groupbit) {
	(*input->vars)["x0"] = Var("x0", s->x0[ip][0]);
	(*input->vars)["y0"] = Var("y0", s->x0[ip][1]);
	(*input->vars)["z0"] = Var("z0", s->x0[ip][2]);
	(*input->vars)["x"] = Var("x", s->x[ip][0]);
	(*input->vars)["y"] = Var("y", s->x[ip][1]);
	(*input->vars)["z"] = Var("z", s->x[ip][2]);
	  
	if (s_set[0])  s->sigma[ip](0,0) = s_value[0].result(mpm);
	if (s_set[1])  s->sigma[ip](1,1) = s_value[1].result(mpm);
	if (s_set[2])  s->sigma[ip](2,2) = s_value[2].result(mpm);
	if (s_set[3])  s->sigma[ip](1,2) = s->sigma[ip](2,1) = s_value[3].result(mpm);
	if (s_set[4])  s->sigma[ip](0,2) = s->sigma[ip](2,0) = s_value[4].result(mpm);
	if (s_set[5])  s->sigma[ip](0,1) = s->sigma[ip](1,0) = s_value[5].result(mpm);

	if (tl) {
	  s->vol0PK1[ip] = s->vol0[ip] * s->sigma[ip];
	}
      }
    }
  }
}
