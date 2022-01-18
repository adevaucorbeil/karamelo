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

#include "fix_neumann_bc_mech.h"
#include "domain.h"
#include "error.h"
#include "group.h"
#include "input.h"
#include "special_functions.h"
#include "universe.h"
#include "update.h"
#include <Eigen/Eigen>
#include <iostream>
#include <vector>

using namespace std;
using namespace FixConst;
using namespace Eigen;


FixNeumannBCMech::FixNeumannBCMech(MPM *mpm, vector<string> args) : Fix(mpm, args)
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
    error->all(FLERR, "Error: too few arguments for fix_heat_flux_particles.\n" +
                          usage.find(domain->dimension)->second);
  }

  if (group->pon[igroup].compare("particles") !=0 ) {
    error->one(FLERR, "fix_convection_nodes needs to be given a group of particles" +
                          group->pon[igroup] + ", " + args[2] +
                          " is a group of " + group->pon[igroup] + ".\n");
  }
  if (universe->me == 0) {
    cout << "Creating new fix FixNeumannBCMech with ID: " << args[0] << endl;
  }

  id = args[0];

  t[0] = input->parsev(args[3]);

  if (domain->dimension >= 2) {
    t[1] = input->parsev(args[4]);
  }
  if (domain->dimension >= 3) {
    t[2] = input->parsev(args[5]);
  }
}

FixNeumannBCMech::~FixNeumannBCMech()
{
}

void FixNeumannBCMech::init()
{
}

void FixNeumannBCMech::setup()
{
}

void FixNeumannBCMech::setmask() {
  mask = 0;
  mask |= INITIAL_INTEGRATE;
}


void FixNeumannBCMech::initial_integrate() {
  // Go through all the particles in the group and set v_update to the right value:
  int solid = group->solid[igroup];
  Solid *s;

  int n = 0;
  double Ap = 0;

  Eigen::Vector3d f, ftot, ftot_reduced;

  f.setZero();
  ftot.setZero();

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

	  if (domain->dimension == 1)
	    Ap = 1;
	  else if (domain->dimension == 2)
	    Ap = sqrt(s->vol[ip]);
	  else 		
	    Ap = pow(s->vol[ip], 2/3);


	  for (int i = 0; i < domain->dimension; i++) {
            f[i] = Ap * t[i].result(mpm);
          }
          s->mbp[ip] += f;
          ftot += f;
        }
        n++;
      }
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

	if (domain->dimension == 1)
	  Ap = 1;
	else if (domain->dimension == 2)
	  Ap = sqrt(s->vol[ip]);
	else 		
	  Ap = pow(s->vol[ip], 2/3);

	for (int i = 0; i < domain->dimension; i++) {
	  f[i] = Ap * t[i].result(mpm);
	}
	s->mbp[ip] += f;
	ftot += f;
      }
    }
 }

  // Reduce ftot:
  MPI_Allreduce(ftot.data(), ftot_reduced.data(), 3, MPI_DOUBLE, MPI_SUM,
                universe->uworld);

  (*input->vars)[id + "_x"] = Var(id + "_x", ftot_reduced[0]);
  (*input->vars)[id + "_y"] = Var(id + "_y", ftot_reduced[1]);
  (*input->vars)[id + "_z"] = Var(id + "_z", ftot_reduced[2]);
}

void FixNeumannBCMech::write_restart(ofstream *of) {
  for (int i = 0; i < domain->dimension; i++)
    t[i].write_to_restart(of);
}

void FixNeumannBCMech::read_restart(ifstream *ifr) {
  for (int i = 0; i < domain->dimension; i++)
    t[i].read_from_restart(ifr);
}
