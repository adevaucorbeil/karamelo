#include <iostream>
#include <vector>
#include "fix_check_solution.h"
#include "input.h"
#include "group.h"
#include "domain.h"
#include "input.h"
#include "update.h"
#include "output.h"
#include "math_special.h"
#include <Eigen/Eigen>

using namespace std;
using namespace FixConst;
using namespace MathSpecial;
using namespace Eigen;

FixChecksolution::FixChecksolution(MPM *mpm, vector<string> args) : Fix(mpm, args)
{
  if (args.size() < 6) {
    cout << "Error: too few arguments for fix_check_solution: requires at least 6 arguments. " << args.size() << " received" << endl;
    exit(1);
  }

  if (group->pon[igroup].compare("particles") !=0 && group->pon[igroup].compare("all") !=0) {
    cout << "fix_check_solution needs to be given a group of nodes" << group->pon[igroup] << ", " << args[2] << " is a group of "<< group->pon[igroup] << "." << endl;
    exit(1);
  }
  cout << "Creating new fix FixChecksolution with ID: " << args[0] << endl;
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

FixChecksolution::~FixChecksolution()
{
}

void FixChecksolution::init()
{
}

void FixChecksolution::setup()
{
}

void FixChecksolution::setmask() {
  mask = 0;
  mask |= FINAL_INTEGRATE;
}


void FixChecksolution::final_integrate() {
  if (update->ntimestep != output->next && update->ntimestep != update->nsteps) return;
  // cout << "In FixChecksolution::post_particles_to_grid()\n";

  // Go through all the nodes in the group and set b to the right value:
  double ux, uy, uz;
    
  int solid = group->solid[igroup];

  int nmax;
  int *mask;
  double *mass;
  double *vol0;
  Eigen::Vector3d error;
  Eigen::Vector3d *x0;
  Eigen::Vector3d *x;  

  error.setZero();

  double vtot;

  if (solid == -1) {
    vtot = 0;
    for (int isolid = 0; isolid < domain->solids.size(); isolid++) {
      vtot += domain->solids[isolid]->vtot;
      x0 = domain->solids[isolid]->x0;
      x = domain->solids[isolid]->x;
      vol0 = domain->solids[isolid]->vol0;
      nmax = domain->solids[isolid]->np;
      mask = domain->solids[isolid]->mask;

      for (int in = 0; in < nmax; in++) {
	if (mask[in] & groupbit) {
	  (*input->vars)["x"] = Var("x", x0[in][0]);
	  (*input->vars)["y"] = Var("y", x0[in][1]);
	  (*input->vars)["z"] = Var("z", x0[in][2]);
	  if (xset) {
	    ux = xvalue.result(mpm);
	    error[0] += vol0[in]*square(ux-(x[in][0]-x0[in][0]));
	  }
	  if (yset) {
	    uy = yvalue.result(mpm);
	    error[1] += vol0[in]*square(uy-(x[in][1]-x0[in][1]));
	  }
	  if (zset) {
	    uz = zvalue.result(mpm);
	    error[2] += vol0[in]*square(uz-(x[in][2]-x0[in][2]));
	  }
	}
      }
    }

    (*input->vars)[id]=Var(id, sqrt((error[0] + error[1] + error[2])/vtot));
    // cout << "f for " << n << " nodes from solid " << domain->solids[isolid]->id << " set." << endl;
  } else {
    vtot = domain->solids[solid]->vtot;
    x0 = domain->solids[solid]->x0;
    x = domain->solids[solid]->x;
    vol0 = domain->solids[solid]->vol0;
    nmax = domain->solids[solid]->np;
    mask = domain->solids[solid]->mask;

    for (int in = 0; in < nmax; in++) {
      if (mask[in] & groupbit) {
	(*input->vars)["x"] = Var("x", x0[in][0]);
	(*input->vars)["y"] = Var("y", x0[in][1]);
	(*input->vars)["z"] = Var("z", x0[in][2]);
	if (xset) {
	  ux = xvalue.result(mpm);
	  error[0] += vol0[in]*square(ux-(x[in][0]-x0[in][0]));
	}
	if (yset) {
	  uy = yvalue.result(mpm);
	  error[1] += vol0[in]*square(uy-(x[in][1]-x0[in][1]));
	}
	if (zset) {
	  uz = zvalue.result(mpm);
	  error[2] += vol0[in]*square(uz-(x[in][2]-x0[in][2]));
	}
      }
    }

    (*input->vars)[id]=Var(id, sqrt((error[0] + error[1] + error[2])/vtot));
    // cout << "f for " << n << " nodes from solid " << domain->solids[solid]->id << " set." << endl;
  }
  // cout << "ftot = [" << ftot[0] << ", " << ftot[1] << ", " << ftot[2] << "]\n"; 
}
