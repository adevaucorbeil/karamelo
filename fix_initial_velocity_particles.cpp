#include <iostream>
#include <vector>
#include "fix_initial_velocity_particles.h"
#include "input.h"
#include "group.h"
#include "domain.h"
#include "update.h"
#include <Eigen/Eigen>

using namespace std;
using namespace FixConst;
using namespace Eigen;

FixInitialVelocityParticles::FixInitialVelocityParticles(MPM *mpm, vector<string> args) : Fix(mpm, args)
{
  if (args.size() < 6) {
    cout << "Error: too few arguments for fix_initial_velocity_particles: requires at least 6 arguments. " << args.size() << " received" << endl;
    exit(1);
  }

  if (group->pon[igroup].compare("particles") !=0 && group->pon[igroup].compare("all") !=0) {
    cout << "fix_initial_velocity_particles needs to be given a group of particles" << group->pon[igroup] << ", " << args[2] << " is a group of "<< group->pon[igroup] << "." << endl;
    exit(1);
  }
  cout << "Creating new fix FixInitialVelocityParticles with ID: " << args[0] << endl;
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

FixInitialVelocityParticles::~FixInitialVelocityParticles()
{
}

void FixInitialVelocityParticles::init()
{
}

void FixInitialVelocityParticles::setup()
{
}

void FixInitialVelocityParticles::setmask() {
  mask = 0;
  mask |= INITIAL_INTEGRATE;
}


void FixInitialVelocityParticles::initial_integrate() {
  if (update->ntimestep !=1) return;
  // cout << "In FixInitialVelocityParticles::initial_integrate()" << endl;

  // Go through all the particles in the group and set v to the right value:
  double vx, vy, vz;
  
  int solid = group->solid[igroup];

  Eigen::Vector3d *v;
  int nmax;
  int *mask;
  int n = 0;
  Eigen::Vector3d *x;

  if (solid == -1) {
    for (int isolid = 0; isolid < domain->solids.size(); isolid++) {
      n = 0;
      v = domain->solids[isolid]->v;
      x = domain->solids[isolid]->x;
      nmax = domain->solids[isolid]->np;
      mask = domain->solids[isolid]->mask;

      for (int ip = 0; ip < nmax; ip++) {
	if (mask[ip] & groupbit) {
	  (*input->vars)["x"] = Var("x", x[ip][0]);
	  (*input->vars)["y"] = Var("y", x[ip][1]);
	  (*input->vars)["z"] = Var("z", x[ip][2]);
	  if (xset) {
	    vx = xvalue.result(mpm);
	    v[ip][0] = vx;
	  }
	  if (yset) {
	    vy = yvalue.result(mpm);
	    v[ip][1] = vy;
	  }
	  if (zset) {
	    vz = zvalue.result(mpm);
	    v[ip][2] = vz;
	  }
	  n++;
	}
      }
      // cout << "v_update for " << n << " particles from solid " << domain->solids[isolid]->id << " set." << endl;
    }
  } else {

    v = domain->solids[solid]->v;
    x = domain->solids[solid]->x;
    nmax = domain->solids[solid]->np;
    mask = domain->solids[solid]->mask;

    for (int ip = 0; ip < nmax; ip++) {
      if (mask[ip] & groupbit) {
	(*input->vars)["x"] = Var("x", x[ip][0]);
	(*input->vars)["y"] = Var("y", x[ip][1]);
	(*input->vars)["z"] = Var("z", x[ip][2]);
	if (xset) {
	  vx = xvalue.result(mpm);
	  v[ip][0] = vx;
	}
	if (yset) {
	  vy = yvalue.result(mpm);
	  v[ip][1] = vy;
	}
	if (zset) {
	  vz = zvalue.result(mpm);
	  v[ip][2] = vz;
	}
	n++;
      }
    }
    // cout << "v_update for " << n << " particles from solid " << domain->solids[solid]->id << " set." << endl;
  }
}
