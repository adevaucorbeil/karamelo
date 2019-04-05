#include <iostream>
#include <vector>
#include "fix_force_nodes.h"
#include "input.h"
#include "group.h"
#include "domain.h"
#include <Eigen/Eigen>

using namespace std;
using namespace FixConst;
using namespace Eigen;

FixForceNodes::FixForceNodes(MPM *mpm, vector<string> args) : Fix(mpm, args)
{
  if (args.size() < 6) {
    cout << "Error: too few arguments for fix_force_nodes: requires at least 6 arguments. " << args.size() << " received" << endl;
    exit(1);
  }

  if (group->pon[igroup].compare("nodes") !=0 ) {
    cout << "fix_force_nodes needs to be given a group of nodes" << group->pon[igroup] << ", " << args[2] << " is a group of "<< group->pon[igroup] << "." << endl;
    exit(1);
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

  Eigen::Vector3d *b;
  int nmax;
  int *mask;
  double *mass;
  int n = 0;
  Eigen::Vector3d ftot;
  ftot.setZero();

  if (solid == -1) {
    for (int isolid = 0; isolid < domain->solids.size(); isolid++) {
      n = 0;
      b = domain->solids[isolid]->grid->b;
      nmax = domain->solids[isolid]->grid->nnodes;
      mask = domain->solids[isolid]->grid->mask;
      mass = domain->solids[isolid]->grid->mass;

      for (int in = 0; in < nmax; in++) {
	if (mass[in] > 0) {
	  if (mask[in] & groupbit) {
	    if (xset) {
	      b[in][0] += fx;
	      ftot[0] += fx;
	    }
	    if (yset) {
	      b[in][1] += fy;
	      ftot[1] += fy;
	    }
	    if (zset) {
	      b[in][2] += fz;
	      ftot[2] += fz;
	    }
	    n++;
	  }
	}
      }
      // cout << "f for " << n << " nodes from solid " << domain->solids[isolid]->id << " set." << endl;
    }
  } else {

    b = domain->solids[solid]->grid->b;
    nmax = domain->solids[solid]->grid->nnodes;
    mask = domain->solids[solid]->grid->mask;
    mass = domain->solids[solid]->grid->mass;

    for (int in = 0; in < nmax; in++) {
      if (mass[in] > 0) {
	if (mask[in] & groupbit) {
	  if (xset) {
	    b[in][0] += fx;
	    ftot[0] += fx;
	  }
	  if (yset) {
	    b[in][1] += fy;
	    ftot[1] += fy;
	  }
	  if (zset) {
	    b[in][2] += fz;
	    ftot[2] += fz;
	  }
	  n++;
	}
      }
    }
    // cout << "f for " << n << " nodes from solid " << domain->solids[solid]->id << " set." << endl;
  }
  // cout << "ftot = [" << ftot[0] << ", " << ftot[1] << ", " << ftot[2] << "]\n"; 
}
