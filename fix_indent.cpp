#include <iostream>
#include <vector>
#include <string>
#include "fix_indent.h"
#include "input.h"
#include "group.h"
#include "domain.h"
#include "input.h"
#include <Eigen/Eigen>
#include "error.h"

using namespace std;
using namespace FixConst;
using namespace Eigen;

FixIndent::FixIndent(MPM *mpm, vector<string> args) : Fix(mpm, args)
{
  if (args.size() < 9) {
    error->all(FLERR,"Error: too few arguments for fix_body_force_current_config: requires at least 9 arguments. " + to_string(args.size()) + " received.\n");
  }

  if (group->pon[igroup].compare("particles") !=0 && group->pon[igroup].compare("all") !=0) {
    error->all(FLERR, "fix_indent needs to be given a group of nodes" + group->pon[igroup] + ", " + args[2] + " is a group of " + group->pon[igroup] + ".\n");
  }
  cout << "Creating new fix FixIndent with ID: " << args[0] << endl;
  id = args[0];

  type = args[3];
  if (args[3].compare("sphere")==0) {
    type = "sphere";
  } else {
    error->all(FLERR,"Error indent type " + args[3] + " unknown. Only type sphere is supported.\n");
  }

  Kpos = 4;
  xpos = 5;
  ypos = 6;
  zpos = 7;
  Rpos = 8;
}

FixIndent::~FixIndent()
{
}

void FixIndent::init()
{
}

void FixIndent::setup()
{
}

void FixIndent::setmask() {
  mask = 0;
  mask |= INITIAL_INTEGRATE;
}


void FixIndent::initial_integrate() {
  // cout << "In FixIndent::initial_integrate()\n";

  // Go through all the particles in the group and set b to the right value:
  Eigen::Vector3d f;

  int solid = group->solid[igroup];

  int nmax;
  int *mask;
  double *mass;
  Eigen::Vector3d ftot;
  Eigen::Vector3d *x;
  Eigen::Vector3d *mb;

  double K = input->parsev(args[Kpos]).result(mpm);
  double R = input->parsev(args[Rpos]).result(mpm);
  Eigen::Vector3d xs(input->parsev(args[xpos]).result(mpm),
		     input->parsev(args[ypos]).result(mpm),
		     input->parsev(args[zpos]).result(mpm));
  Eigen::Vector3d xsp;

  double r, dr, fmag;

  ftot.setZero();

  if (solid == -1) {
    for (int isolid = 0; isolid < domain->solids.size(); isolid++) {
      mb = domain->solids[isolid]->mb;
      x = domain->solids[isolid]->x;
      nmax = domain->solids[isolid]->np;
      mask = domain->solids[isolid]->mask;
      mass = domain->solids[isolid]->mass;

      for (int ip = 0; ip < nmax; ip++) {
	if (mass[ip] > 0) {
	  if (mask[ip] & groupbit) {
	    // Gross screening:
	    xsp = x[ip] - xs;

	    if (( xsp[0] < R ) && ( xsp[1] < R ) && ( xsp[2] < R )
		&& ( xsp[0] > -R ) && ( xsp[1] > -R ) && ( xsp[2] > -R )) {

	      r = xsp.norm();
	      // Finer screening:
	      if (r < R) {
		dr = r - R;
		fmag = K*dr*dr;
		// Maybe fmag should be inversely proportional to the mass of the particle!!
		f = fmag*xsp/r;
		mb[ip] += f;
		ftot += f;
	      }
	    }
	  }
	}
      }
      (*input->vars)[id+"_x"]=Var(id+"_x", ftot[0]);
      (*input->vars)[id+"_y"]=Var(id+"_y", ftot[1]);
      (*input->vars)[id+"_z"]=Var(id+"_z", ftot[2]);
      // cout << "f for " << n << " nodes from solid " << domain->solids[isolid]->id << " set." << endl;
    }
  } else {
    mb = domain->solids[solid]->mb;
    x = domain->solids[solid]->x;
    nmax = domain->solids[solid]->np;
    mask = domain->solids[solid]->mask;
    mass = domain->solids[solid]->mass;

    for (int ip = 0; ip < nmax; ip++) {
      if (mass[ip] > 0) {
	if (mask[ip] & groupbit) {
	  // Gross screening:
	  xsp = x[ip] - xs;
	  if (( xsp[0] < R ) && ( xsp[1] < R ) && ( xsp[2] < R )
	      && ( xsp[0] > -R ) && ( xsp[1] > -R ) && ( xsp[2] > -R )) {

	    r = xsp.norm();
	    // Finer screening:
	    if (r < R) {
	      dr = r - R;
	      fmag = K*dr*dr;
	      // Maybe fmag should be inversely proportional to the mass of the particle!!
	      f = fmag*xsp/r;
	      mb[ip] += f;
	      ftot += f;
	    }
	  }
	}
      }
    }

    (*input->vars)[id+"_x"]=Var(id+"_x", ftot[0]);
    (*input->vars)[id+"_y"]=Var(id+"_y", ftot[1]);
    (*input->vars)[id+"_z"]=Var(id+"_z", ftot[2]);
    // cout << "f for " << n << " nodes from solid " << domain->solids[isolid]->id << " set." << endl;
    }
  // cout << "ftot = [" << ftot[0] << ", " << ftot[1] << ", " << ftot[2] << "]\n"; 
}
