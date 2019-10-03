#include <iostream>
#include <vector>
#include <string>
#include <Eigen/Eigen>
#include "fix_velocity_nodes.h"
#include "input.h"
#include "group.h"
#include "domain.h"
#include "grid.h"
#include "error.h"

using namespace std;
using namespace FixConst;
using namespace Eigen;

FixVelocityNodes::FixVelocityNodes(MPM *mpm, vector<string> args) : Fix(mpm, args)
{
  if (domain->dimension == 3 && args.size()<6) {
    error->all(FLERR,"Error: too few arguments for fix_velocity_nodes: requires at least 6 arguments. " + to_string(args.size()) + " received.\n");
  } else if (domain->dimension == 2 && args.size()<5) {
    error->all(FLERR,"Error: too few arguments for fix_velocity_nodes: requires at least 5 arguments. " + to_string(args.size()) + " received.\n");
  } else if (domain->dimension == 1 && args.size()<4) {
    error->all(FLERR,"Error: too few arguments for fix_velocity_nodes: requires at least 4 arguments. " + to_string(args.size()) + " received.\n");
  }

  if (group->pon[igroup].compare("nodes") !=0 ) {
    error->all(FLERR, "fix_velocity_nodes needs to be given a group of nodes" + group->pon[igroup] + ", " + args[2] + " is a group of "+ group->pon[igroup] + ".\n");
  }
  cout << "Creating new fix FixVelocityNodes with ID: " << args[0] << endl;
  id = args[0];

  xset = yset = zset = false;

  args_previous_step = args;

  string time = "time";

  if (args[3].compare("NULL") != 0) {
    // xvalue = input->parsev(args[3]);
    xpos = 3;
    xset = true;

    // Replace "time" by "time - dt" in the x argument:
    while(args_previous_step[xpos].find(time)!=std::string::npos) {
      args_previous_step[xpos].replace(args_previous_step[xpos].find(time),time.length(),"time - dt");
    }
  }

  if (domain->dimension >= 2) {
    if (args[4].compare("NULL") != 0) {
      ypos = 4;
      // yvalue = input->parsev(args[4]);
      yset = true;
      // Replace "time" by "time - dt" in the y argument:
      while(args_previous_step[ypos].find(time)!=std::string::npos) {
	args_previous_step[ypos].replace(args_previous_step[ypos].find(time),time.length(),"time - dt");
      }
    }
  }

  if (domain->dimension == 3) {
    if (args[5].compare("NULL") != 0) {
      zpos = 5;
      // zvalue = input->parsev(args[5]);
      zset = true;

      // Replace "time" by "time - dt" in the z argument:
      while(args_previous_step[zpos].find(time)!=std::string::npos) {
	args_previous_step[zpos].replace(args_previous_step[zpos].find(time),time.length(),"time - dt");
      }
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
  double vx_old, vy_old, vz_old;

  if (xset) {
    vx = input->parsev(args[xpos]).result(mpm);
    vx_old = input->parsev(args_previous_step[xpos]).result(mpm);
    // cout << "Set v_update[0] to " << xvalue.eq() << "=" << vx << endl;
    // cout << "Set v[0] to " << vx_old << endl;
  }

  if (yset) {
    vy = input->parsev(args[ypos]).result(mpm);
    vy_old = input->parsev(args_previous_step[ypos]).result(mpm);
    // cout << "Set v_update[1] to " << "=" <<  vy << endl;
    // cout << "Set v[1] to " << "=" <<  vy_old << endl;
  }

  if (zset) {
    vz = input->parsev(args[zpos]).result(mpm);
    vz_old = input->parsev(args_previous_step[zpos]).result(mpm);
    // cout << "Set v_update[2] to " << "=" <<  vz << endl;
    // cout << "Set v[2] to " << "=" <<  vz_old << endl;
  }

  int solid = group->solid[igroup];
  Grid *g;

  if (solid == -1) {
    for (int isolid = 0; isolid < domain->solids.size(); isolid++) {
      g = domain->solids[isolid]->grid;

      for (int ip = 0; ip < g->nnodes_local + g->nnodes_ghost; ip++) {
	if (g->mask[ip] & groupbit) {
	  if (xset) {
	    g->v_update[ip][0] = vx;
	    g->v[ip][0] = vx_old;
	  }
	  if (yset) {
	    g->v_update[ip][1] = vy;
	    g->v[ip][1] = vy_old;
	  }
	  if (zset) {
	    g->v_update[ip][2] = vz;
	    g->v[ip][2] = vz_old;
	  }
	}
      }
    }
  } else {

    g = domain->solids[solid]->grid;

    for (int ip = 0; ip < g->nnodes_local + g->nnodes_ghost; ip++) {
      if (g->mask[ip] & groupbit) {
	if (xset) {
	  g->v_update[ip][0] = vx;
	  g->v[ip][0] = vx_old;
	}
	if (yset) {
	  g->v_update[ip][1] = vy;
	  g->v[ip][1] = vy_old;
	}
	if (zset) {
	  g->v_update[ip][2] = vz;
	  g->v[ip][2] = vz_old;
	}
      }
    }
  }
}

void FixVelocityNodes::post_velocities_to_grid() {
  // cout << "In FixVelocityNodes::post_velocities_to_grid()" << endl;

  // Go through all the particles in the group and set v to the right value:
  double vx, vy, vz;

  if (xset) {
    vx = input->parsev(args[xpos]).result(mpm);
  }

  if (yset) {
    vy = input->parsev(args[ypos]).result(mpm);
  }

  if (zset) {
    vz = input->parsev(args[zpos]).result(mpm);
  }
  
  int solid = group->solid[igroup];
  Grid *g;

  if (solid == -1) {
    for (int isolid = 0; isolid < domain->solids.size(); isolid++) {
      g = domain->solids[isolid]->grid;

      for (int ip = 0; ip < g->nnodes_local + g->nnodes_ghost; ip++) {
	if (g->mask[ip] & groupbit) {
	  if (xset) g->v[ip][0] = vx;
	  if (yset) g->v[ip][1] = vy;
	  if (zset) g->v[ip][2] = vz;
	}
      }
    }
  } else {
    g = domain->solids[solid]->grid;

    for (int ip = 0; ip < g->nnodes_local + g->nnodes_ghost; ip++) {
      if (g->mask[ip] & groupbit) {
	if (xset) g->v[ip][0] = vx;
	if (yset) g->v[ip][1] = vy;
	if (zset) g->v[ip][2] = vz;
      }
    }
  }
}
