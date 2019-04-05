#include <iostream>
#include <vector>
#include "fix_indent.h"
#include "input.h"
#include "group.h"
#include "domain.h"
#include <Eigen/Eigen>

using namespace std;
using namespace FixConst;
using namespace Eigen;

FixIndent::FixIndent(MPM *mpm, vector<string> args) : Fix(mpm, args)
{
  if (args.size() < 6) {
    cout << "Error: too few arguments for fix_indent: requires at least 8 arguments. " << args.size() << " received" << endl;
    exit(1);
  }

  if (group->pon[igroup].compare("nodes") !=0 && group->pon[igroup].compare("all") !=0) {
    cout << "fix_indent needs to be given a group of nodes" << group->pon[igroup] << ", " << args[2] << " is a group of "<< group->pon[igroup] << "." << endl;
    exit(1);
  }
  cout << "Creating new fix FixIndent with ID: " << args[0] << endl;
  id = args[0];

  xvalue = input->parsev(args[3]);
  yvalue = input->parsev(args[4]);
  zvalue = input->parsev(args[5]);
  R = input->parsev(args[6]);
  K = input->parsev(args[7]);
}

void FixIndent::setmask() {
  mask = 0;
  mask |= POST_PARTICLES_TO_GRID;
}


void FixIndent::post_particles_to_grid() {
  // cout << "In FixIndent::post_particles_to_grid()\n";

  // Go through all the nodes in the group and set b to the right value:
  Eigen::Vector3d xind;

  xind[0] = xvalue.result(mpm);
  xind[1] = yvalue.result(mpm);
  xind[2] = zvalue.result(mpm);

    
  int solid = group->solid[igroup];

  Eigen::Vector3d *b;
  Eigen::Vector3d *xn;
  Eigen::Vector3d f;
  Eigen::Vector3d delx;
  int nmax;
  int *mask;
  double *mass;
  int n = 0;
  double r, dr, fmag;

  if (solid == -1) {
    for (int isolid = 0; isolid < domain->solids.size(); isolid++) {
      n = 0;
      b = domain->solids[isolid]->grid->b;
      xn = domain->solids[isolid]->grid->x;
      nmax = domain->solids[isolid]->grid->nnodes;
      mask = domain->solids[isolid]->grid->mask;
      mass = domain->solids[isolid]->grid->mass;

      for (int in = 0; in < nmax; in++) {
	if (mass[in] > 0) {
	  if (mask[in] & groupbit) {
	    delx = xn[in] - xind;
	    r = delx.norm();
	    dr = R - r;

	    if (dr > 0.0) {
	      fmag = K*dr*dr;
	      f = delx*fmag/r;
	      b[in] += f;
	      // cout << "b[" << in <<"]=[" << b[in][0] << ", " << b[in][1] << ", " << b[in][2] << "]\n";
	      domain->solids[isolid]->dtCFL = MIN(mass[in] / (dr * K), domain->solids[isolid]->dtCFL);
	    }
	    n++;
	  }
	}
      }
      // cout << "b for " << n << " nodes from solid " << domain->solids[isolid]->id << " set." << endl;
    }
  } else {

    b = domain->solids[solid]->grid->b;
    xn = domain->solids[solid]->grid->x;
    nmax = domain->solids[solid]->grid->nnodes;
    mask = domain->solids[solid]->grid->mask;
    mass = domain->solids[solid]->grid->mass;

    for (int in = 0; in < nmax; in++) {
      if (mass[in] > 0) {
	if (mask[in] & groupbit) {
	  delx = xn[in] - xind;
	  r = delx.norm();
	  dr = R - r;

	  if (dr > 0.0) {
	    fmag = K*dr*dr;
	    f = delx*fmag/r;
	    b[in] += f;
	    // cout << "b[" << in <<"]=[" << b[in][0] << ", " << b[in][1] << ", " << b[in][2] << "]\n";
	    domain->solids[solid]->dtCFL = MIN(mass[in] / (dr * K), domain->solids[solid]->dtCFL);
	  }
	  n++;
	}
      }
    }
    // cout << "b for " << n << " nodes from solid " << domain->solids[solid]->id << " set." << endl;
  }
}
