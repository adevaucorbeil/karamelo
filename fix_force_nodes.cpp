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
#include <Eigen/Eigen>
#include "fix_force_nodes.h"
#include "input.h"
#include "group.h"
#include "domain.h"
#include "input.h"
#include "universe.h"
#include "grid.h"
#include "error.h"

using namespace std;
using namespace FixConst;
using namespace Eigen;

FixForceNodes::FixForceNodes(MPM *mpm, vector<string> args) : Fix(mpm, args)
{
  if (args.size() < 6) {
    error->all(FLERR,"Error: too few arguments for fix_force_nodes: requires at least 6 arguments. " + to_string(args.size()) + " received.\n");
  }

  if (group->pon[igroup].compare("nodes") !=0 ) {
    error->all(FLERR,"fix_force_nodes needs to be given a group of nodes" + group->pon[igroup] + ", " + args[2] + " is a group of " + group->pon[igroup] + ".\n");
  }
  cout << "Creating new fix FixForceNodes with ID: " << args[0] << endl;
  id = args[0];

  xset = yset = zset = false;

  if (args[3].compare("NULL") != 0) {
    xvalue = input->parsev(args[3]);
    xset = true;
  }

  if (args[4].compare("NULL") != 0) {
    yvalue = input->parsev(args[4]);
    yset = true;
  }

  if (args[5].compare("NULL") != 0) {
    zvalue = input->parsev(args[5]);
    zset = true;
  }
}

FixForceNodes::~FixForceNodes()
{
}

void FixForceNodes::init()
{
}

void FixForceNodes::setup()
{
}

void FixForceNodes::setmask() {
  mask = 0;
  mask |= POST_PARTICLES_TO_GRID;
}


void FixForceNodes::post_particles_to_grid() {
  // cout << "In FixForceNodes::post_particles_to_grid()\n";

  // Go through all the nodes in the group and set b to the right value:
  double fx, fy, fz;

  if (xset) {
    fx = xvalue.result(mpm);
    // cout << "Set v_update[0] to " << xvalue.eq() << "=" << fx << endl;
  }

  if (yset) {
    fy = yvalue.result(mpm);
    // cout << "Set v_update[1] to " << "=" <<  fy << endl;
  }

  if (zset) {
    fz = zvalue.result(mpm);
    // cout << "Set v_update[2] to " << "=" <<  fz << endl;
  }
    
  int solid = group->solid[igroup];
  Grid *g;

  int n = 0;
  Eigen::Vector3d ftot, ftot_reduced;
  ftot.setZero();

  if (solid == -1) {
    for (int isolid = 0; isolid < domain->solids.size(); isolid++) {
      g = domain->solids[isolid]->grid;
      n = 0;

      for (int in = 0; in < g->nnodes_local + g->nnodes_ghost; in++) {
	if (g->mass[in] > 0) {
	  if (g->mask[in] & groupbit) {
	    n++;
	  }
	}
      }
	    
      for (int in = 0; in < g->nnodes_local + g->nnodes_ghost; in++) {
	if (g->mass[in] > 0) {
	  if (g->mask[in] & groupbit) {
	    if (xset) {
	      g->mb[in][0] += fx/((double) n);
	      if(in < g->nnodes_local) ftot[0] += fx/((double) n);
	    }
	    if (yset) {
	      g->mb[in][1] += fy/((double) n);
	      if(in < g->nnodes_local) ftot[1] += fy/((double) n);
	    }
	    if (zset) {
	      g->mb[in][2] += fz/((double) n);
	      if(in < g->nnodes_local) ftot[2] += fz/((double) n);
	    }
	  }
	}
      }
    }
  } else {

    g = domain->solids[solid]->grid;
    n = 0;

    for (int in = 0; in < g->nnodes_local + g->nnodes_ghost; in++) {
      if (g->mass[in] > 0) {
	if (g->mask[in] & groupbit) {
	  n++;
	}
      }
    }

    for (int in = 0; in < g->nnodes_local + g->nnodes_ghost; in++) {
      if (g->mass[in] > 0) {
	if (g->mask[in] & groupbit) {
	  if (xset) {
	    g->mb[in][0] += fx/((double) n);
	    if(in < g->nnodes_local) ftot[0] += fx/((double) n);
	  }
	  if (yset) {
	    g->mb[in][1] += fy/((double) n);
	    if(in < g->nnodes_local) ftot[1] += fy/((double) n);
	  }
	  if (zset) {
	    g->mb[in][2] += fz/((double) n);
	    if(in < g->nnodes_local) ftot[2] += fz/((double) n);
	  }
	}
      }
    }
  }

  // Reduce ftot:
  MPI_Allreduce(ftot.data(),ftot_reduced.data(),3,MPI_DOUBLE,MPI_SUM,universe->uworld);

  if (xset) (*input->vars)[id+"_x"]=Var(id+"_x", ftot_reduced[0]);
  if (yset) (*input->vars)[id+"_y"]=Var(id+"_y", ftot_reduced[1]);
  if (zset) (*input->vars)[id+"_z"]=Var(id+"_z", ftot_reduced[2]);
  // cout << "f for " << n << " nodes from solid " << domain->solids[solid]->id << " set." << endl;
  // cout << "ftot = [" << ftot[0] << ", " << ftot[1] << ", " << ftot[2] << "]\n"; 
}
