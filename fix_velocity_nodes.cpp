#include <iostream>
#include <vector>
#include "fix_velocity_nodes.h"
#include "input.h"
#include "group.h"
#include "domain.h"
#include <Eigen/Eigen>

using namespace std;
using namespace FixConst;
using namespace Eigen;

FixVelocityNodes::FixVelocityNodes(MPM *mpm, vector<string> args) : Fix(mpm, args)
{
  if (domain->dimension == 3 && args.size()<6) {
    cout << "Error: too few arguments for fix_velocity_nodes: requires at least 6 arguments. " << args.size() << " received" << endl;
    exit(1);
  } else if (domain->dimension == 2 && args.size()<5) {
    cout << "Error: too few arguments for fix_velocity_nodes: requires at least 5 arguments. " << args.size() << " received" << endl;
    exit(1);
  } else if (domain->dimension == 1 && args.size()<4) {
    cout << "Error: too few arguments for fix_velocity_nodes: requires at least 4 arguments. " << args.size() << " received" << endl;
    exit(1);
  }

  if (group->pon[igroup].compare("nodes") !=0 ) {
    cout << "fix_velocity_nodes needs to be given a group of nodes" << group->pon[igroup] << ", " << args[2] << " is a group of "<< group->pon[igroup] << "." << endl;
    exit(1);
  }
  cout << "Creating new fix FixVelocityNodes with ID: " << args[0] << endl;
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

FixVelocityNodes::~FixVelocityNodes()
{
}

void FixVelocityNodes::init()
{
}

void FixVelocityNodes::setup()
{
}

void FixVelocityNodes::setmask() {
  mask = 0;
  mask |= POST_UPDATE_GRID_STATE;
  mask |= POST_VELOCITIES_TO_GRID;
}


void FixVelocityNodes::post_update_grid_state() {
  // cout << "In FixVelocityNodes::post_update_grid_state()" << endl;

  // Go through all the nodes in the group and set v_update to the right value:
  double vx, vy, vz;

  if (xset) {
    vx = xvalue.result(mpm);
    // cout << "Set v_update[0] to " << xvalue.eq() << "=" << vx << endl;
  }

  if (yset) {
    vy = yvalue.result(mpm);
    // cout << "Set v_update[1] to " << "=" <<  vy << endl;
  }

  if (zset) {
    vz = zvalue.result(mpm);
    // cout << "Set v_update[2] to " << "=" <<  vz << endl;
  }

  
  int solid = group->solid[igroup];

  Eigen::Vector3d *v_update;
  int nmax;
  int *mask;
  int n = 0;

  if (solid == -1) {
    for (int isolid = 0; isolid < domain->solids.size(); isolid++) {
      n = 0;
      v_update = domain->solids[isolid]->grid->v_update;
      nmax = domain->solids[isolid]->grid->nnodes;
      mask = domain->solids[isolid]->grid->mask;

      for (int ip = 0; ip < nmax; ip++) {
	if (mask[ip] & groupbit) {
	  if (xset) v_update[ip][0] = vx;
	  if (yset) v_update[ip][1] = vy;
	  if (zset) v_update[ip][2] = vz;
	  n++;
	}
      }
      // cout << "v_update for " << n << " nodes from solid " << domain->solids[isolid]->id << " set." << endl;
    }
  } else {

    v_update = domain->solids[solid]->grid->v_update;
    nmax = domain->solids[solid]->grid->nnodes;
    mask = domain->solids[solid]->grid->mask;

    for (int ip = 0; ip < nmax; ip++) {
      if (mask[ip] & groupbit) {
	if (xset) v_update[ip][0] = vx;
	if (yset) v_update[ip][1] = vy;
	if (zset) v_update[ip][2] = vz;
	n++;
      }
    }
    // cout << "v_update for " << n << " nodes from solid " << domain->solids[solid]->id << " set." << endl;
  }
}

void FixVelocityNodes::post_velocities_to_grid() {
  // cout << "In FixVelocityNodes::post_velocities_to_grid()" << endl;

  // Go through all the particles in the group and set v to the right value:
  double vx, vy, vz;

  if (xset) {
    vx = xvalue.result(mpm);
    // cout << "Set v[0] to " << xvalue.eq() << "=" << vx << endl;
  }

  if (yset) {
    vy = yvalue.result(mpm);
    // cout << "Set v[1] to " << "=" <<  vy << endl;
  }

  if (zset) {
    vz = zvalue.result(mpm);
    // cout << "Set v[2] to " << "=" <<  vz << endl;
  }

  
  int solid = group->solid[igroup];

  Eigen::Vector3d *v;
  int nmax;
  int *mask;
  int n = 0;

  if (solid == -1) {
    for (int isolid = 0; isolid < domain->solids.size(); isolid++) {
      n = 0;
      v = domain->solids[isolid]->grid->v;
      nmax = domain->solids[isolid]->grid->nnodes;
      mask = domain->solids[isolid]->grid->mask;

      for (int ip = 0; ip < nmax; ip++) {
	if (mask[ip] & groupbit) {
	  if (xset) v[ip][0] = vx;
	  if (yset) v[ip][1] = vy;
	  if (zset) v[ip][2] = vz;
	  n++;
	}
      }
      // cout << "v for " << n << " nodes from solid " << domain->solids[isolid]->id << " set." << endl;
    }
  } else {

    v = domain->solids[solid]->grid->v;
    nmax = domain->solids[solid]->grid->nnodes;
    mask = domain->solids[solid]->grid->mask;

    for (int ip = 0; ip < nmax; ip++) {
      if (mask[ip] & groupbit) {
	if (xset) v[ip][0] = vx;
	if (yset) v[ip][1] = vy;
	if (zset) v[ip][2] = vz;
	n++;
      }
    }
    // cout << "v for " << n << " nodes from solid " << domain->solids[solid]->id << " set." << endl;
  }
}
