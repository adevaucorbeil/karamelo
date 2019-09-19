#include <iostream>
#include <vector>
#include "fix_body_force_current_config.h"
#include "input.h"
#include "group.h"
#include "domain.h"
#include "input.h"
#include <Eigen/Eigen>

using namespace std;
using namespace FixConst;
using namespace Eigen;

FixBodyforceCurrentConfig::FixBodyforceCurrentConfig(MPM *mpm, vector<string> args) : Fix(mpm, args)
{
  if (args.size() < 6) {
    cout << "Error: too few arguments for fix_body_force_current_config: requires at least 6 arguments. " << args.size() << " received" << endl;
    exit(1);
  }

  if (group->pon[igroup].compare("nodes") !=0 && group->pon[igroup].compare("all") !=0) {
    cout << "fix_body_force_current_config needs to be given a group of nodes" << group->pon[igroup] << ", " << args[2] << " is a group of "<< group->pon[igroup] << "." << endl;
    exit(1);
  }
  cout << "Creating new fix FixBodyforceCurrentConfig with ID: " << args[0] << endl;
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

FixBodyforceCurrentConfig::~FixBodyforceCurrentConfig()
{
}

void FixBodyforceCurrentConfig::init()
{
}

void FixBodyforceCurrentConfig::setup()
{
}

void FixBodyforceCurrentConfig::setmask() {
  mask = 0;
  mask |= INITIAL_INTEGRATE;
}


void FixBodyforceCurrentConfig::initial_integrate() {
  // cout << "In FixBodyforceCurrentConfig::initial_integrate()\n";

  // Go through all the nodes in the group and set b to the right value:
  Eigen::Vector3d f;
    
  int solid = group->solid[igroup];

  Eigen::Vector3d *mb;
  int nmax;
  int *mask;
  double *mass;
  Eigen::Vector3d ftot;
  Eigen::Vector3d *x0;
  Eigen::Matrix3d *R;

  ftot.setZero();

  if (solid == -1) {
    for (int isolid = 0; isolid < domain->solids.size(); isolid++) {
      mb = domain->solids[isolid]->mb;
      x0 = domain->solids[isolid]->x0;
      nmax = domain->solids[isolid]->np;
      mask = domain->solids[isolid]->mask;
      mass = domain->solids[isolid]->mass;
      R = domain->solids[isolid]->R;

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
	  }
	}
      }
      if (xset) (*input->vars)[id+"_x"]=Var(id+"_x", ftot[0]);
      if (yset) (*input->vars)[id+"_y"]=Var(id+"_y", ftot[1]);
      if (zset) (*input->vars)[id+"_z"]=Var(id+"_z", ftot[2]);
      // cout << "f for " << n << " nodes from solid " << domain->solids[isolid]->id << " set." << endl;
    }
  } else {

    mb = domain->solids[solid]->mb;
    x0 = domain->solids[solid]->x0;
    nmax = domain->solids[solid]->np;
    mask = domain->solids[solid]->mask;
    mass = domain->solids[solid]->mass;
    R = domain->solids[solid]->R;

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
	}
      }
    }
    if (xset) (*input->vars)[id+"_x"]=Var(id+"_x", ftot[0]);
    if (yset) (*input->vars)[id+"_y"]=Var(id+"_y", ftot[1]);
    if (zset) (*input->vars)[id+"_z"]=Var(id+"_z", ftot[2]);
    // cout << "f for " << n << " nodes from solid " << domain->solids[isolid]->id << " set." << endl;
    }
  // cout << "ftot = [" << ftot[0] << ", " << ftot[1] << ", " << ftot[2] << "]\n"; 
}
