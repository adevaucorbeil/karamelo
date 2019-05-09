#include <iostream>
#include <vector>
#include "fix_body_force.h"
#include "input.h"
#include "group.h"
#include "domain.h"
#include "input.h"
#include <Eigen/Eigen>

using namespace std;
using namespace FixConst;
using namespace Eigen;

FixBodyforce::FixBodyforce(MPM *mpm, vector<string> args) : Fix(mpm, args)
{
  if (args.size() < 6) {
    cout << "Error: too few arguments for fix_body_force: requires at least 6 arguments. " << args.size() << " received" << endl;
    exit(1);
  }

  if (group->pon[igroup].compare("nodes") !=0 && group->pon[igroup].compare("all") !=0) {
    cout << "fix_body_force needs to be given a group of nodes" << group->pon[igroup] << ", " << args[2] << " is a group of "<< group->pon[igroup] << "." << endl;
    exit(1);
  }
  cout << "Creating new fix FixBodyforce with ID: " << args[0] << endl;
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
  double fx, fy, fz;
    
  int solid = group->solid[igroup];

  Eigen::Vector3d *b;
  int nmax;
  int *mask;
  double *mass;
  Eigen::Vector3d ftot;
  Eigen::Vector3d *x0;

  ftot.setZero();

  if (solid == -1) {
    for (int isolid = 0; isolid < domain->solids.size(); isolid++) {
      b = domain->solids[isolid]->grid->b;
      x0 = domain->solids[isolid]->grid->x0;
      nmax = domain->solids[isolid]->grid->nnodes;
      mask = domain->solids[isolid]->grid->mask;
      mass = domain->solids[isolid]->grid->mass;

      for (int in = 0; in < nmax; in++) {
	if (mass[in] > 0) {
	  if (mask[in] & groupbit) {
	      (*input->vars)["x"] = Var("x", x0[in][0]);
	      (*input->vars)["y"] = Var("y", x0[in][1]);
	      (*input->vars)["z"] = Var("z", x0[in][2]);
	    if (xset) {
	      fx = xvalue.result(mpm);
	      b[in][0] += mass[in]*fx;
	      ftot[0] += mass[in]*fx;
	    }
	    if (yset) {
	      fy = yvalue.result(mpm);
	      b[in][1] += mass[in]*fy;
	      ftot[1] += mass[in]*fy;
	    }
	    if (zset) {
	      fz = zvalue.result(mpm);
	      b[in][2] += mass[in]*fz;
	      ftot[2] += mass[in]*fz;
	    }
	  }
	}
      }
      if (xset) (*input->vars)[id+"_x"]=Var(id+"_x", ftot[0]);
      if (yset) (*input->vars)[id+"_y"]=Var(id+"_y", ftot[1]);
      if (zset) (*input->vars)[id+"_z"]=Var(id+"_z", ftot[2]);
      // cout << "f for " << n << " nodes from solid " << domain->solids[isolid]->id << " set." << endl;
    }
  } else {

    b = domain->solids[solid]->grid->b;
    x0 = domain->solids[solid]->grid->x0;
    nmax = domain->solids[solid]->grid->nnodes;
    mask = domain->solids[solid]->grid->mask;
    mass = domain->solids[solid]->grid->mass;

    
    for (int in = 0; in < nmax; in++) {
      if (mass[in] > 0) {
	if (mask[in] & groupbit) {
	  (*input->vars)["x"] = Var("x", x0[in][0]);
	  (*input->vars)["y"] = Var("y", x0[in][1]);
	  (*input->vars)["z"] = Var("z", x0[in][2]);
	  if (xset) {
	    fx = xvalue.result(mpm);
	    b[in][0] += mass[in]*fx;
	    ftot[0] += mass[in]*fx;
	  }
	  if (yset) {
	    fy = yvalue.result(mpm);
	    b[in][1] += mass[in]*fy;
	    ftot[1] += mass[in]*fy;
	  }
	  if (zset) {
	    fz = zvalue.result(mpm);
	    b[in][2] += mass[in]*fz;
	    ftot[2] += mass[in]*fz;
	  }
	}
      }
    }
    if (xset) (*input->vars)[id+"_x"]=Var(id+"_x", ftot[0]);
    if (yset) (*input->vars)[id+"_y"]=Var(id+"_y", ftot[1]);
    if (zset) (*input->vars)[id+"_z"]=Var(id+"_z", ftot[2]);
    // cout << "f for " << n << " nodes from solid " << domain->solids[solid]->id << " set." << endl;
  }
  // cout << "ftot = [" << ftot[0] << ", " << ftot[1] << ", " << ftot[2] << "]\n"; 
}
