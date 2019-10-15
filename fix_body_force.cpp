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
  if (domain->dimension == 3 && args.size()<6) {
    cout << "Error: too few arguments for fix_body_force: requires at least 6 arguments. " << args.size() << " received" << endl;
    exit(1);
  } else if (domain->dimension == 2 && args.size()<5) {
    cout << "Error: too few arguments for fix_body_force: requires at least 5 arguments. " << args.size() << " received" << endl;
    exit(1);
  } else if (domain->dimension == 1 && args.size()<4) {
    cout << "Error: too few arguments for fix_body_force: requires at least 4 arguments. " << args.size() << " received" << endl;
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
  Eigen::Vector3d f;

  int solid = group->solid[igroup];

  Eigen::Vector3d *mb;
  int nmax;
  int *mask;
  double *mass;
  Eigen::Vector3d ftot;
  Eigen::Vector3d *x0;
  Eigen::Vector3d *x;
  Eigen::Matrix3d *R;

  // double mtot = 0;
  ftot.setZero();

  if (solid == -1) {
    for (int isolid = 0; isolid < domain->solids.size(); isolid++) {
      mb = domain->solids[isolid]->grid->mb;
      x0 = domain->solids[isolid]->grid->x0;
      x = domain->solids[isolid]->grid->x;
      nmax = domain->solids[isolid]->grid->nnodes;
      mask = domain->solids[isolid]->grid->mask;
      mass = domain->solids[isolid]->grid->mass;

      for (int in = 0; in < nmax; in++) {
	if (mass[in] > 0) {
	  if (mask[in] & groupbit) {
	      (*input->vars)["x"] = Var("x", x[in][0]);
	      (*input->vars)["y"] = Var("y", x[in][1]);
	      (*input->vars)["z"] = Var("z", x[in][2]);
	      (*input->vars)["x0"] = Var("x0", x[in][0]);
	      (*input->vars)["y0"] = Var("y0", x[in][1]);
	      (*input->vars)["z0"] = Var("z0", x[in][2]);

	      f.setZero();
	      if (xset) f[0] = xvalue.result(mpm);
	      if (yset) f[1] = yvalue.result(mpm);
	      if (zset) f[2] = zvalue.result(mpm);

	      f *= mass[in];
	      mb[in] += f;
	      ftot += f;
	      // mtot += mass[in];
	  }
	}
      }
      if (xset) (*input->vars)[id+"_x"]=Var(id+"_x", ftot[0]);
      if (yset) (*input->vars)[id+"_y"]=Var(id+"_y", ftot[1]);
      if (zset) (*input->vars)[id+"_z"]=Var(id+"_z", ftot[2]);
      // cout << "f for " << n << " nodes from solid " << domain->solids[isolid]->id << " set." << endl;
    }
  } else {

    mb = domain->solids[solid]->grid->mb;
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

	  f.setZero();
	  if (xset) f[0] = xvalue.result(mpm);
	  if (yset) f[1] = yvalue.result(mpm);
	  if (zset) f[2] = zvalue.result(mpm);

	  f *= mass[in];
	  mb[in] += f;
	  ftot += f;
	  // mtot += mass[in];
	}
      }
    }
    if (xset) (*input->vars)[id+"_x"]=Var(id+"_x", ftot[0]);
    if (yset) (*input->vars)[id+"_y"]=Var(id+"_y", ftot[1]);
    if (zset) (*input->vars)[id+"_z"]=Var(id+"_z", ftot[2]);
    // cout << "f for " << n << " nodes from solid " << domain->solids[solid]->id << " set." << endl;
  }
  // cout << "ftot = [" << ftot[0] << ", " << ftot[1] << ", " << ftot[2] << "], mass = " << mtot << "\n"; 
}
