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

#include <iostream>
#include <vector>
#include <string>
#include <matrix.h>
#include <fix_body_force.h>
#include <input.h>
#include <group.h>
#include <domain.h>
#include <input.h>
#include <universe.h>
#include <grid.h>
#include <error.h>

using namespace std;
using namespace FixConst;


FixBodyforce::FixBodyforce(MPM *mpm, vector<string> args) : Fix(mpm, args)
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

    xset = yset = zset = false;
    return;
  }

  if (domain->dimension == 3 && args.size()<6) {
    error->all(FLERR,"Error: too few arguments for fix_body_force: requires at least 6 arguments. " + to_string(args.size()) + " received.\n");
  } else if (domain->dimension == 2 && args.size()<5) {
    error->all(FLERR,"Error: too few arguments for fix_body_force: requires at least 5 arguments. " + to_string(args.size()) + " received.\n");
  } else if (domain->dimension == 1 && args.size()<4) {
    error->all(FLERR,"Error: too few arguments for fix_body_force: requires at least 4 arguments. " + to_string(args.size()) + " received.\n");
  }

  if (group->pon[igroup].compare("nodes") !=0 && group->pon[igroup].compare("all") !=0) {
    error->all(FLERR,"fix_body_force needs to be given a group of nodes" + group->pon[igroup] + ", " + args[2] + " is a group of " + group->pon[igroup] + ".\n");
  }
  if (universe->me == 0) {
    cout << "Creating new fix FixBodyforce with ID: " << args[0] << endl;
  }
  id = args[0];

  xset = yset = zset = false;

  if (args[3].compare("NULL") != 0) {
    xvalue = input->parsev(args[3]);
    xset = true;
  }

  if (domain->dimension >= 2) {
    if (args[4].compare("NULL") != 0) {
      yvalue = input->parsev(args[4]);
      yset = true;
    }
  }

  if (domain->dimension == 3) {
    if (args[5].compare("NULL") != 0) {
      zvalue = input->parsev(args[5]);
      zset = true;
    }
  }
}

FixBodyforce::~FixBodyforce()
{
}

void FixBodyforce::init()
{
}

void FixBodyforce::setup()
{
}

void FixBodyforce::setmask() {
  mask = 0;
  mask |= POST_PARTICLES_TO_GRID;
}


void FixBodyforce::post_particles_to_grid() {
  // cout << "In FixBodyforce::post_particles_to_grid()\n";

  // Go through all the nodes in the group and set b to the right value:
  Vector3d f;

  int solid = group->solid[igroup];
  Grid *g;

  Vector3d ftot, ftot_reduced;

  // double mtot = 0;
  ftot = Vector3d();

  if (solid == -1) {
    for (int isolid = 0; isolid < domain->solids.size(); isolid++) {
      g = domain->solids[isolid]->grid;

      for (int in = 0; in < g->nnodes_local + g->nnodes_ghost; in++) {
	if (g->mass[in] > 0) {
	  if (g->mask[in] & groupbit) {
	      (*input->vars)["x0"] = Var("x0", g->x0[in][0]);
	      (*input->vars)["y0"] = Var("y0", g->x0[in][1]);
	      (*input->vars)["z0"] = Var("z0", g->x0[in][2]);
              f = Vector3d();
	      if (xset) f[0] = xvalue.result(mpm);
	      if (yset) f[1] = yvalue.result(mpm);
	      if (zset) f[2] = zvalue.result(mpm);

	      f *= g->mass[in];
	      g->mb[in] += f;
	      if (in < g->nnodes_local) {
		ftot += f;
		// mtot += g->mass[in];
	      }
	  }
	}
      }
    }
  } else {
    g = domain->solids[solid]->grid;

    for (int in = 0; in < g->nnodes_local + g->nnodes_ghost; in++) {
      if (g->mass[in] > 0) {
	if (g->mask[in] & groupbit) {
	  (*input->vars)["x0"] = Var("x0", g->x0[in][0]);
	  (*input->vars)["y0"] = Var("y0", g->x0[in][1]);
	  (*input->vars)["z0"] = Var("z0", g->x0[in][2]);

	  f = Vector3d();
	  if (xset) f[0] = xvalue.result(mpm);
	  if (yset) f[1] = yvalue.result(mpm);
	  if (zset) f[2] = zvalue.result(mpm);

	  f *= g->mass[in];
	  g->mb[in] += f;

	  if (in < g->nnodes_local) {
	    ftot += f;
	    // mtot += g->mass[in];
	  }
	}
      }
    }
  }

  // Reduce ftot:
  MPI_Allreduce(ftot.elements,ftot_reduced.elements,3,MPI_DOUBLE,MPI_SUM,universe->uworld);

  if (xset) (*input->vars)[id+"_x"]=Var(id+"_x", ftot_reduced[0]);
  if (yset) (*input->vars)[id+"_y"]=Var(id+"_y", ftot_reduced[1]);
  if (zset) (*input->vars)[id+"_z"]=Var(id+"_z", ftot_reduced[2]);
  // cout << "ftot = [" << ftot[0] << ", " << ftot[1] << ", " << ftot[2] << "], mass = " << mtot << "\n"; 
}

void FixBodyforce::write_restart(ofstream *of) {
  of->write(reinterpret_cast<const char *>(&xset), sizeof(bool));
  of->write(reinterpret_cast<const char *>(&yset), sizeof(bool));
  of->write(reinterpret_cast<const char *>(&zset), sizeof(bool));

  if (xset)
    xvalue.write_to_restart(of);
  if (yset)
    yvalue.write_to_restart(of);
  if (zset)
    zvalue.write_to_restart(of);
}

void FixBodyforce::read_restart(ifstream *ifr) {
  ifr->read(reinterpret_cast<char *>(&xset), sizeof(bool));
  ifr->read(reinterpret_cast<char *>(&yset), sizeof(bool));
  ifr->read(reinterpret_cast<char *>(&zset), sizeof(bool));

  if (xset)
    xvalue.read_from_restart(ifr);
  if (yset)
    yvalue.read_from_restart(ifr);
  if (zset)
    zvalue.read_from_restart(ifr);
}
