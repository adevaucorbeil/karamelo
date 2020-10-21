/* ----------------------------------------------------------------------
 *
 *                    ***       Karamelo       ***
 *               Parallel Material Point Method Simulator
 *
 * Copyright (2019) Alban de Vaucorbeil, alban.devaucorbeil@monash.edu
 * Materials Science and Engineering, Monash University
 * Clayton VIC 3800, Australia

 * This software is distributed under the GNU General Public License.
 *
 * ----------------------------------------------------------------------- */

#include "fix_meld_tool.h"
#include "domain.h"
#include "error.h"
#include "group.h"
#include "input.h"
#include "solid.h"
#include "universe.h"
#include "update.h"
#include <Eigen/Eigen>
#include <iostream>
#include <string>
#include <vector>

using namespace std;
using namespace FixConst;
using namespace Eigen;

FixMeldTool::FixMeldTool(MPM *mpm, vector<string> args)
    : Fix(mpm, args) {
  if (args.size() < Nargs) {
    error->all(FLERR, "Error: not enough arguments.\n" + usage);
  }

  if (group->pon[igroup].compare("particles") != 0 &&
      group->pon[igroup].compare("all") != 0) {
    error->all(FLERR, "fix_meldtool needs to be given a group of nodes" +
                          group->pon[igroup] + ", " + args[2] +
                          " is a group of " + group->pon[igroup] + ".\n");
  }
  cout << "Creating new fix FixMeldTool with ID: " << args[0] << endl;
  id = args[0];
  K = input->parsev(args[3]).result(mpm);
  if (args[4].compare("x") == 0) {
    dim = X;
    axis0 = Y;
    axis1 = Z;
  } else if (args[4].compare("y") == 0) {
    dim = Y;
    axis0 = X;
    axis1 = Z;
  } else if (args[4].compare("z") == 0) {
    dim = Z;
    axis0 = X;
    axis1 = Y;
  } else {
    error->all(FLERR, "Unknown dim: " + args[4] + "! Options are : x, y, and z.\n");
  }

  w = input->parsev(args[5]).result(mpm);
  lo = input->parsev(args[9]).result(mpm);
  hi = input->parsev(args[10]).result(mpm);
  Rmax = input->parsev(args[11]).result(mpm);
  RmaxSq = Rmax * Rmax;
}

FixMeldTool::~FixMeldTool() {}

void FixMeldTool::init() {}

void FixMeldTool::setup() {}

void FixMeldTool::setmask() {
  mask = 0;
  mask |= INITIAL_INTEGRATE;
}

void FixMeldTool::initial_integrate() {
  // cout << "In FixMeldTool::initial_integrate()\n";

  // Go through all the particles in the group and set b to the right value:
  Eigen::Vector3d f;

  int solid = group->solid[igroup];

  Solid *s;
  Eigen::Vector3d ftot, ftot_reduced;


  double theta = input->parsev(args[thetapos]).result(mpm);
  double c1 = input->parsev(args[c1pos]).result(mpm);
  double c2 = input->parsev(args[c2pos]).result(mpm);

  Eigen::Vector3d xprime, c;
  Eigen::Matrix3d R, Rt;

  if (dim == X) {
    R(0,0) = 1;
    R(0,1) = 0;
    R(0,2) = 0;
    R(1,0) = 0;
    R(1,1) = cos(theta);
    R(1,2) = sin(theta);
    R(2,0) = 0;
    R(2,1) = -sin(theta);
    R(2,2) = cos(theta);

    c(0) = lo;
    c(1) = c1;
    c(2) = c2;
  } else if (dim == Y) {
    R(0,0) = cos(theta);
    R(0,1) = 0;
    R(0,2) = sin(theta);
    R(1,0) = 0;
    R(1,1) = 1;
    R(1,2) = 0;
    R(2,0) = -sin(theta);
    R(2,1) = 0;
    R(2,2) = cos(theta);

    c(0) = c1;
    c(1) = lo;
    c(2) = c2;
  } else if (dim == Z) {
    R(0,0) = 1;
    R(0,1) = 0;
    R(0,2) = 0;
    R(1,0) = 0;
    R(1,1) = cos(theta);
    R(1,2) = sin(theta);
    R(2,0) = 0;
    R(2,1) = -sin(theta);
    R(2,2) = cos(theta);

    c(0) = c1;
    c(1) = c2;
    c(2) = lo;
  }

  Rt = R.transpose();

  ftot.setZero();

  double p0, p1, p2, p, rSq, fmag;

  if (solid == -1) {
    for (int isolid = 0; isolid < domain->solids.size(); isolid++) {
      s = domain->solids[isolid];
      double KG = K * s->mat->G;

      for (int ip = 0; ip < s->np_local; ip++) {
	if (s->mass[ip] > 0) {
	  if (s->mask[ip] & groupbit) {

	    if (s->x[ip][dim] >= lo && s->x[ip][dim] <= hi &&
		s->x[ip](axis0) <= Rmax && s->x[ip](axis0) >= -Rmax &&
		s->x[ip](axis1) <= Rmax && s->x[ip](axis1) >= -Rmax) {

	      rSq = s->x[ip](axis0) * s->x[ip](axis0) + s->x[ip](axis1) * s->x[ip](axis1);

	      if (rSq <= RmaxSq) {
		xprime = R * (s->x[ip] - c);
		p0 = xprime[axis0];
		p1 = xprime[axis1];
		p2 = xprime[dim];

		f.setZero();

		if (p0 > w) {
		  p = p0 - w;
		  f[axis0] = -p;
		}

		if (p0 < -w) {
		  p = -w - p0;
		  f[axis0] = p;
		}

		if (p1 > w) {
		  p = p1 - w;
		  f[axis1] = -p;
		}

		if (p1 < -w) {
		  p = -p1 - w;
		  f[axis1] = p;
		}

		if (p2 > 0) {
		  if (p2 < fabs(f[axis0]) || p2 < fabs(f[axis1])) {
		    f.setZero();
		    f[dim] = p2;
		  }
		}

		f = KG * (1.0 - s->damage[ip]) * (Rt * f);
		s->mbp[ip] += f;
		ftot += f;
	      }
	    }
          }
        }
      }
    }
  } else {
    s = domain->solids[solid];
    double KG = K * s->mat->G;

    for (int ip = 0; ip < s->np_local; ip++) {
      if (s->mass[ip] > 0) {
	if (s->mask[ip] & groupbit) {

	  if (s->x[ip][dim] >= lo && s->x[ip][dim] <= hi &&
	      s->x[ip](axis0) <= Rmax && s->x[ip](axis0) >= -Rmax &&
	      s->x[ip](axis1) <= Rmax && s->x[ip](axis1) >= -Rmax) {

	    rSq = s->x[ip](axis0) * s->x[ip](axis0) + s->x[ip](axis1) * s->x[ip](axis1);

	    if (rSq <= RmaxSq) {
	      xprime = R * (s->x[ip] - c);
	      p0 = xprime[axis0];
	      p1 = xprime[axis1];
	      p2 = xprime[dim];

	      f.setZero();

	      if (p0 > w) {
		p = p0 - w;
		f[axis0] = -p;
	      }

	      if (p0 < -w) {
		p = -w - p0;
		f[axis0] = p;
	      }

	      if (p1 > w) {
		p = p1 - w;
		f[axis1] = -p;
	      }

	      if (p1 < -w) {
		p = -p1 - w;
		f[axis1] = p;
	      }

	      if (p2 > 0) {
		if (p2 < fabs(f[axis0]) || p2 < fabs(f[axis1])) {
		  f.setZero();
		  f[dim] = p2;
		}
	      }

	      f = KG * (1.0 - s->damage[ip]) * (Rt * f);
	      s->mbp[ip] += f;
	      ftot += f;
	    }
	  }
	}
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
