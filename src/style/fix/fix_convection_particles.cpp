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

#include <fix_convection_particles.h>
#include <domain.h>
#include <error.h>
#include <group.h>
#include <input.h>
#include <special_functions.h>
#include <universe.h>
#include <update.h>
#include <iostream>
#include <vector>

using namespace std;
using namespace FixConst;


FixConvectionParticles::FixConvectionParticles(MPM *mpm, vector<string> args) : Fix(mpm, args)
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
    error->all(FLERR, "Error: too few arguments for fix_convection_particles.\n" +
                          usage);
  }

  if (args.size() > Nargs) {
    error->all(FLERR, "Error: many few arguments for fix_convection_particles.\n" +
                          usage);
  }

  if (group->pon[igroup].compare("particles") !=0 ) {
    error->one(FLERR, "fix_convection_nodes needs to be given a group of particles" +
                          group->pon[igroup] + ", " + args[2] +
                          " is a group of " + group->pon[igroup] + ".\n");
  }
  if (universe->me == 0) {
    cout << "Creating new fix FixConvectionParticles with ID: " << args[0] << endl;
  }

  id = args[0];

  h = input->parsev(args[3]);
  Tinf = input->parsev(args[4]);
}

FixConvectionParticles::~FixConvectionParticles()
{
}

void FixConvectionParticles::init()
{
}

void FixConvectionParticles::setup()
{
}

void FixConvectionParticles::setmask() {
  mask = 0;
  mask |= INITIAL_INTEGRATE;
  mask |= FINAL_INTEGRATE;
}


void FixConvectionParticles::initial_integrate() {
  // Go through all the particles in the group and set v_update to the right value:
  int solid = group->solid[igroup];
  Solid *s;

  double qtemp, q, qtot, qtot_reduced, invcp, Ap = 0;

  if (solid == -1) {
    for (int isolid = 0; isolid < domain->solids.size(); isolid++) {
      s = domain->solids[isolid];
      invcp = s->mat->invcp;

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

	  qtemp = h * (Tinf.result(mpm) - s->T[ip]);
	  s->gamma[ip] += Ap * qtemp * invcp;
	  qtot += qtemp;
	}
      }
    }
  } else {
    s = domain->solids[solid];
    invcp = s->mat->invcp;

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

	qtemp = h * (Tinf.result(mpm) - s->T[ip]);
	s->gamma[ip] += Ap * qtemp * invcp;
	qtot += qtemp;
      }
    }
  }

  // Reduce qtot:
  MPI_Allreduce(&qtot, &qtot_reduced, 1, MPI_DOUBLE, MPI_SUM,
                universe->uworld);

  (*input->vars)[id + "_s"] = Var(id + "_s", qtot_reduced);
}

void FixConvectionParticles::write_restart(ofstream *of) {
  of->write(reinterpret_cast<const char *>(&h), sizeof(double));
  Tinf.write_to_restart(of);
}

void FixConvectionParticles::read_restart(ifstream *ifr) {
  ifr->read(reinterpret_cast<char *>(&h), sizeof(double));
  Tinf.read_from_restart(ifr);
}
