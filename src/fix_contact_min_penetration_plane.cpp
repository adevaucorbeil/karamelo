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

#include "fix_contact_min_penetration_plane.h"
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

#define four_thirds 1.333333333

FixContactMinPenetrationPlane::FixContactMinPenetrationPlane(MPM *mpm, vector<string> args)
    : Fix(mpm, args) {
  if (args.size() < Nargs) {
    error->all(FLERR, "Error: not enough arguments.\n" + usage);
  }

  solid = domain->find_solid(args[2]);
  if (solid < 0) {
    error->all(FLERR, "Error: solid " + args[2] + " unknown.\n");
  }

  cout << "Creating new fix FixContactMinPenetrationPlane with ID: " << args[0] << endl;
  id = args[0];
  requires_ghost_particles = false;

  xq[0] = input->parsev(args[3]);
  xq[1] = input->parsev(args[4]);
  xq[2] = input->parsev(args[5]);

  n[0] = input->parsev(args[6]);
  n[1] = input->parsev(args[7]);
  n[2] = input->parsev(args[8]);

  mu = input->parsev(args[9]);

  // Normalize:
  n /= n.norm();
  
  D = -n[0] * xq[0] - n[1] * xq[1] - n[2] * xq[2];
}

FixContactMinPenetrationPlane::~FixContactMinPenetrationPlane() {}

void FixContactMinPenetrationPlane::init() {}

void FixContactMinPenetrationPlane::setup() {}

void FixContactMinPenetrationPlane::setmask() {
  mask = 0;
  mask |= INITIAL_INTEGRATE;
  //mask |= POST_ADVANCE_PARTICLES;
}
void FixContactMinPenetrationPlane::initial_integrate() {
  // cout << "In FixContactMinPenetrationPlane::initial_integrate()\n";

  // Go through all the particles in the group and set b to the right value:
  Eigen::Vector3d f;

  Solid *s;
  Eigen::Vector3d ftot, ftot_reduced, vt;

  double Rp, r, p, d, fnorm, vtnorm;

  ftot.setZero();

  s = domain->solids[solid];

  if (domain->dimension == 2) {
    for (int ip = 0; ip < s->np_local; ip++) {
      d = n.dot(s->x[ip] - xq);

      // Extremely gross screening:
      if (d < s->grid->cellsize) {
	Rp = 0.5 * sqrt(s->vol[ip]);

	p = Rp - d;
	// Fine screening:
	if (p >= 0) {
	  fnorm = s->mass[ip] * p / (update->dt * update->dt);
	  f = fnorm * n;

	  if (mu != 0) {
	    vt = s->v[ip] - n.dot(s->v[ip]) * n;
	    vtnorm = vt.norm();
	    if (vtnorm != 0) {
	      vt /= vtnorm;
	      f -= mu * fnorm * vt;
	    }
	  }
	  s->mbp[ip] += f;
	  ftot += f;
	}
      }
    }
  } else if (domain->dimension == 3) {
    for (int ip = 0; ip < s->np_local; ip++) {
      d = n.dot(s->x[ip] - xq);

      // Extremely gross screening:
      if (d < s->grid->cellsize) {
	Rp = 0.5 * cbrt(s->vol[ip]);

	p = Rp - d;
	// Fine screening:
	if (p >= 0) {
	  fnorm = s->mass[ip] * p / (update->dt * update->dt);
	  f = fnorm * n;

	  if (mu != 0) {
	    vt = s->v[ip] - n.dot(s->v[ip]) * n;
	    vtnorm = vt.norm();
	    if (vtnorm != 0) {
	      vt /= vtnorm;
	      f -= mu * fnorm * vt;
	    }
	  }
	  s->mbp[ip] += f;
	  ftot += f;
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
