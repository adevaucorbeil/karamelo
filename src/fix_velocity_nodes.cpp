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
#include "update.h"
#include "universe.h"
using namespace std;
using namespace FixConst;
using namespace Eigen;

FixVelocityNodes::FixVelocityNodes(MPM *mpm, vector<string> args) : Fix(mpm, args)
{
  if (args.size() < 3) {
    error->all(FLERR, "Error: not enough arguments.\n");
  }

  if (args[2].compare("restart") ==
      0) { // If the keyword restart, we are expecting to have read_restart()
           // launched right after.
    igroup = stoi(args[3]);
    if (igroup == -1) {
      cout << "Could not find group number " << args[3] << endl;
    }
    groupbit = group->bitmask[igroup];
    
    xset = yset = zset = false;
    return;
  }

  if (args.size() < Nargs.find(domain->dimension)->second) {
    error->all(FLERR, "Error: too few arguments for fix_velocity_nodes.\n" +
                          usage.find(domain->dimension)->second);
  }

  if (group->pon[igroup].compare("nodes") !=0 ) {
    error->one(FLERR, "fix_velocity_nodes needs to be given a group of nodes" + group->pon[igroup] + ", " + args[2] + " is a group of "+ group->pon[igroup] + ".\n");
  }
  cout << "Creating new fix FixVelocityNodes with ID: " << args[0] << endl;
  id = args[0];

  xset = yset = zset = false;

  string time = "time";

  if (args[3].compare("NULL") != 0) {
    // xvalue = input->parsev(args[3]);
    xset = true;
    xvalue = input->parsev(args[3]);

    string previous = args[3];

    // Replace "time" by "time - dt" in the x argument:
    while(previous.find(time)!=std::string::npos) {
      previous.replace(previous.find(time),time.length(),"time - dt");
    }
    xprevvalue = input->parsev(previous);
  }

  if (domain->dimension >= 2) {
    if (args[4].compare("NULL") != 0) {
      yvalue = input->parsev(args[4]);
      yset = true;

      string previous = args[4];

      // Replace "time" by "time - dt" in the y argument:
      while(previous.find(time)!=std::string::npos) {
	previous.replace(previous.find(time),time.length(),"time - dt");
      }
      yprevvalue = input->parsev(previous);
    }
  }

  if (domain->dimension == 3) {
    if (args[5].compare("NULL") != 0) {
      zvalue = input->parsev(args[5]);
      zset = true;

      string previous = args[5];

      // Replace "time" by "time - dt" in the z argument:
      while(previous.find(time)!=std::string::npos) {
	previous.replace(previous.find(time),time.length(),"time - dt");
      }
      zprevvalue = input->parsev(previous);
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
    vx = xvalue.result(mpm);
    vx_old = xprevvalue.result(mpm);
    // cout << "Set v_update[0] to " << xvalue.eq() << "=" << vx << endl;
    // cout << "Set v[0] to " << vx_old << endl;
  }

  if (yset) {
    vy = yvalue.result(mpm);
    vy_old = yprevvalue.result(mpm);
    // cout << "Set v_update[1] to " << "=" <<  vy << endl;
    // cout << "Set v[1] to " << "=" <<  vy_old << endl;
  }

  if (zset) {
    vz = zvalue.result(mpm);
    vz_old = zprevvalue.result(mpm);
    // cout << "Set v_update[2] to " << "=" <<  vz << endl;
    // cout << "Set v[2] to " << "=" <<  vz_old << endl;
  }

  int solid = group->solid[igroup];
  Grid *g;

  Eigen::Vector3d Dv, ftot, ftot_reduced;
  ftot.setZero();
  double inv_dt = 1.0/update->dt;

  if (solid == -1) {
    for (int isolid = 0; isolid < domain->solids.size(); isolid++) {
      g = domain->solids[isolid]->grid;

      for (int ip = 0; ip < g->nnodes_local + g->nnodes_ghost; ip++) {
	if (g->mask[ip] & groupbit) {
	  Dv.setZero();
	  if (xset) {
	    Dv[0] = vx - g->v_update[ip][0];
	    g->v_update[ip][0] = vx;
	    g->v[ip][0] = vx_old;
	  }
	  if (yset) {
	    Dv[1] = vy - g->v_update[ip][1];
	    g->v_update[ip][1] = vy;
	    g->v[ip][1] = vy_old;
	  }
	  if (zset) {
	    Dv[2] = vz - g->v_update[ip][2];
	    g->v_update[ip][2] = vz;
	    g->v[ip][2] = vz_old;
	  }
          ftot += (inv_dt * g->mass[ip]) * Dv;
	}
      }
    }
  } else {

    g = domain->solids[solid]->grid;

    for (int ip = 0; ip < g->nnodes_local + g->nnodes_ghost; ip++) {
      if (g->mask[ip] & groupbit) {
	Dv.setZero();
	if (xset) {
	  Dv[0] = vx - g->v_update[ip][0];
	  g->v_update[ip][0] = vx;
	  g->v[ip][0] = vx_old;
	}
	if (yset) {
	  Dv[1] = vy - g->v_update[ip][1];
	  g->v_update[ip][1] = vy;
	  g->v[ip][1] = vy_old;
	}
	if (zset) {
	  Dv[2] = vz - g->v_update[ip][2];
	  g->v_update[ip][2] = vz;
	  g->v[ip][2] = vz_old;
	}
	ftot += (inv_dt * g->mass[ip]) * Dv;
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

void FixVelocityNodes::post_velocities_to_grid() {
  // cout << "In FixVelocityNodes::post_velocities_to_grid()" << endl;

  // Go through all the particles in the group and set v to the right value:
  double vx, vy, vz;

  if (xset) {
    vx = xvalue.result(mpm);
  }

  if (yset) {
    vy = yvalue.result(mpm);
  }

  if (zset) {
    vz = zvalue.result(mpm);
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

void FixVelocityNodes::write_restart(ofstream *of) {
  of->write(reinterpret_cast<const char *>(&xset), sizeof(bool));
  of->write(reinterpret_cast<const char *>(&yset), sizeof(bool));
  of->write(reinterpret_cast<const char *>(&zset), sizeof(bool));

  if (xset) {
    xvalue.write_to_restart(of);
    xprevvalue.write_to_restart(of);
  }
  if (yset) {
    yvalue.write_to_restart(of);
    yprevvalue.write_to_restart(of);
  }
  if (zset) {
    zvalue.write_to_restart(of);
    zprevvalue.write_to_restart(of);
  }
}

void FixVelocityNodes::read_restart(ifstream *ifr) {
  ifr->read(reinterpret_cast<char *>(&xset), sizeof(bool));
   ifr->read(reinterpret_cast<char *>(&yset), sizeof(bool));
  ifr->read(reinterpret_cast<char *>(&zset), sizeof(bool));

  if (xset) {
    xvalue.read_from_restart(ifr);
    xprevvalue.read_from_restart(ifr);
  }
  if (yset) {
    yvalue.read_from_restart(ifr);
    yprevvalue.read_from_restart(ifr);
  }
  if (zset) {
    zvalue.read_from_restart(ifr);
    zprevvalue.read_from_restart(ifr);
  }
}
