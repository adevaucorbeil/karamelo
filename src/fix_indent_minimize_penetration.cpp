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

#include "fix_indent_minimize_penetration.h"
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

FixIndentMinimizePenetration::FixIndentMinimizePenetration(MPM *mpm, vector<string> args)
    : Fix(mpm, args) {
  if (args.size() < Nargs) {
    error->all(FLERR, "Error: not enough arguments.\n" + usage);
  }
  if (args.size() > Nargs) {
    error->all(FLERR, "Error: too many arguments.\n" + usage);
  }

  type_pos = 3;
  Rpos = 4;
  xpos = 5;
  ypos = 6;
  zpos = 7;
  vxpos = 8;
  vypos = 9;
  vzpos = 10;

  if (group->pon[igroup].compare("particles") != 0 &&
      group->pon[igroup].compare("all") != 0) {
    error->all(FLERR, "fix_indent_hertz needs to be given a group of particles" +
                          group->pon[igroup] + ", " + args[2] +
                          " is a group of " + group->pon[igroup] + ".\n");
  }
  cout << "Creating new fix FixIndentMinimizePenetration with ID: " << args[0] << endl;
  id = args[0];

  type = args[type_pos];
  if (args[type_pos].compare("sphere") == 0) {
    type = "sphere";
  } else {
    error->all(FLERR, "Error indent type " + args[type_pos] +
                          " unknown. Only type sphere is supported.\n");
  }
  mu = input->parsev(args[11]);
}

FixIndentMinimizePenetration::~FixIndentMinimizePenetration() {}

void FixIndentMinimizePenetration::init() {}

void FixIndentMinimizePenetration::setup() {}

void FixIndentMinimizePenetration::setmask() {
  mask = 0;
  mask |= INITIAL_INTEGRATE;
}

void FixIndentMinimizePenetration::initial_integrate() {
  // Go through all the particles in the group and set b to the right value:
  Eigen::Vector3d f;

  int solid = group->solid[igroup];

  Solid *s;
  Eigen::Vector3d ftot, ftot_reduced, ffric;

  double R = input->parsev(args[Rpos]).result(mpm);
  Eigen::Vector3d xs(input->parsev(args[xpos]).result(mpm),
                     input->parsev(args[ypos]).result(mpm),
                     input->parsev(args[zpos]).result(mpm));
  Eigen::Vector3d vs(input->parsev(args[vxpos]).result(mpm),
                     input->parsev(args[vypos]).result(mpm),
                     input->parsev(args[vzpos]).result(mpm));
  Eigen::Vector3d xsp, vps, vt;

  double Rs, Rp, r, p, fmag, vtnorm, vndotxsp;

  int n, n_reduced;
  ftot.setZero();
  n = 0;

  if (solid == -1) {
    for (int isolid = 0; isolid < domain->solids.size(); isolid++) {
      s = domain->solids[isolid];

      if (domain->dimension == 2) {
        for (int ip = 0; ip < s->np_local; ip++) {
          if (s->mass[ip] > 0) {
            if (s->mask[ip] & groupbit) {

              // Gross screening:
              xsp = s->x[ip] - xs;

	      if (domain->axisymmetric) {
		Rp = 0.5*sqrt(s->vol0[ip] / s->x0[ip][0]);
	      } else {
		Rp = 0.5*sqrt(s->vol0[ip]);
	      }
              Rs = R + Rp;

              if ((xsp[0] < Rs) && (xsp[1] < Rs) && (xsp[0] > -Rs) &&
                  (xsp[1] > -Rs)) {

                r = xsp.norm();
                // Finer screening:
                if (r < Rs) {
                  p = Rs - r; // penetration

                  if (p > 0) {
		    n++;
		    xsp /= r;
                    fmag = s->mass[ip] * p / (update->dt * update->dt);
                    f = fmag * xsp;// / r;

		    if (mu != 0) {
		      vps = vs - s->v[ip];
		      vt = vps - vps.dot(xsp) * xsp;// / (r*r);
		      vtnorm = vt.norm();
		      if (vtnorm != 0) {
			vt /= vtnorm;
			f += MIN(s->mass[ip] * vtnorm / update->dt, mu * fmag) * vt;
		      }
		    }

		    s->mbp[ip] += f;
                    ftot += f;
                  } else {
                    fmag = 0;
                    f.setZero();
                  }
                }
              }
            }
          }
        }
      }

      if (domain->dimension == 3) {
        for (int ip = 0; ip < s->np_local; ip++) {
          if (s->mass[ip] > 0) {
            if (s->mask[ip] & groupbit) {

             // Gross screening:
              xsp = s->x[ip] - xs;

              Rp = 0.5*cbrt(s->vol0[ip]);
              Rs = R + Rp;

              if ((xsp[0] < Rs) && (xsp[1] < Rs) && (xsp[2] < Rs) &&
                  (xsp[0] > -Rs) && (xsp[1] > -Rs) && (xsp[2] > -Rs)) {

                r = xsp.norm();
                // Finer screening:
                if (r < Rs) {
                  p = Rs - r; // penetration

                  if (p > 0) {
		    n++;
		    xsp /= r;
		    vps = vs - s->v[ip];
                    fmag = s->mass[ip] * p / (update->dt * update->dt);
                    f = fmag * xsp;// / r;
		    if (mu != 0) {
		      vndotxsp = vps.dot(xsp);
		      vt = vps - vndotxsp * xsp;// / (r*r);
		      vtnorm = vt.norm();
		      if (vtnorm != 0) {
			// Reference: DOI 10.1002/nag.2233
			if (vtnorm <= mu * vndotxsp) {
			  ffric = s->mass[ip] * vt/ update->dt;
			  //f += MIN(s->mass[ip] * vtnorm / update->dt, mu * fmag) * vt;
			} else {
			  // vt /= vtnorm;
			  // Either:
			  // ffric = mu * fmag * vt;
			  // Or:
			  ffric = mu * s->mass[ip] * abs(vndotxsp) * vt / (update->dt * vtnorm);
			}
			cout << "Particle " << s->ptag[ip]
			     << " f_friction=[" << ffric[0] <<"," << ffric[1] <<"," << ffric[2] << "] "
			     << "f=[" << f[0] <<"," << f[1] <<"," << f[2] << "] "
			     << "vps=[" << vps[0] <<"," << vps[1] <<"," << vps[2] << "] "
			      << "vt=[" << vt[0] <<"," << vt[1] <<"," << vt[2] << "] "
			//      << "vs=[" << vs[0] <<"," << vs[1] <<"," << vs[2] << "] "
			//      << "vp=[" << s->v[ip][0] <<"," << s->v[ip][1] <<"," << s->v[ip][2] << "] "
			// //      << "ap=[" << s->a[ip][0] <<"," << s->a[ip][1] <<"," << s->a[ip][2] << "] "
			     << "xsp=[" << xsp[0] <<"," << xsp[1] <<"," << xsp[2] << "] "
			     << "dt=" << update->dt << endl;
			// " mass=" << s->mass[ip] << " vtnorm=" << vtnorm << "\n";
			// //   cout << endl;

			f += ffric;
		      }
		    }

		    s->mbp[ip] += f;
                    ftot += f;
                  } else {
                    fmag = 0;
                    f.setZero();
                  }
                }
              }
            }
          }
        }
      }
    }
  } else {
    s = domain->solids[solid];
    if (domain->dimension == 2) {
      for (int ip = 0; ip < s->np_local; ip++) {
	if (s->mass[ip] > 0) {
	  if (s->mask[ip] & groupbit) {

	    // Gross screening:
	    xsp = s->x[ip] - xs;

	    if (domain->axisymmetric) {
	      Rp = 0.5*sqrt(s->vol0[ip] / s->x0[ip][0]);
	    } else {
	      Rp = 0.5*sqrt(s->vol0[ip]);
	    }
	    Rs = R + Rp;

	    if ((xsp[0] < Rs) && (xsp[1] < Rs) && (xsp[0] > -Rs) &&
		(xsp[1] > -Rs)) {

	      r = xsp.norm();
	      // Finer screening:
	      if (r < Rs) {
		p = Rs - r; // penetration

		if (p > 0) {
		  n++;
		  xsp /= r;
		  fmag = s->mass[ip] * p / (update->dt * update->dt);
		  f = fmag * xsp;// / r;

		  if (mu != 0) {
		    vps = vs - s->v[ip];
		    vt = vps - vps.dot(xsp) * xsp;// / (r*r);
		    vtnorm = vt.norm();
		    if (vtnorm != 0) {
		      vt /= vtnorm;
		      f += MIN(s->mass[ip] * vtnorm / update->dt, mu * fmag) * vt;
		    }
		  }

		  s->mbp[ip] += f;
		  ftot += f;
		} else {
		  fmag = 0;
		  f.setZero();
		}
	      }
	    }
	  }
	}
      }
    }

    if (domain->dimension == 3) {
      for (int ip = 0; ip < s->np_local; ip++) {
	if (s->mass[ip] > 0) {
	  if (s->mask[ip] & groupbit) {

	    // Gross screening:
	    xsp = s->x[ip] - xs;

	    Rp = 0.5*cbrt(s->vol0[ip]);
	    Rs = R + Rp;

	    if ((xsp[0] < Rs) && (xsp[1] < Rs) && (xsp[2] < Rs) &&
		(xsp[0] > -Rs) && (xsp[1] > -Rs) && (xsp[2] > -Rs)) {

	      r = xsp.norm();
	      // Finer screening:
	      if (r < Rs) {
		p = Rs - r; // penetration

		if (p > 0) {
		  n++;
		  xsp /= r;
		  fmag = s->mass[ip] * p / (update->dt * update->dt);
		  f = fmag * xsp;// / r;

		  if (mu != 0) {
		    vps = vs - s->v[ip];
		    vt = vps - vps.dot(xsp) * xsp;// / (r*r);
		    vtnorm = vt.norm();
		    if (vtnorm != 0) {
		      vt /= vtnorm;
		      f += MIN(s->mass[ip] * vtnorm / update->dt, mu * fmag) * vt;
		    }
		  }

		  s->mbp[ip] += f;
		  ftot += f;
		} else {
		  fmag = 0;
		  f.setZero();
		}
	      }
	    }
	  }
	}
      }
    }
  }

  // Reduce ftot:
  MPI_Allreduce(&n, &n_reduced, 1, MPI_INT, MPI_SUM, universe->uworld);
  MPI_Allreduce(ftot.data(), ftot_reduced.data(), 3, MPI_DOUBLE, MPI_SUM,
                universe->uworld);

  (*input->vars)[id + "_s"] = Var(id + "_s", n_reduced);
  (*input->vars)[id + "_x"] = Var(id + "_x", ftot_reduced[0]);
  (*input->vars)[id + "_y"] = Var(id + "_y", ftot_reduced[1]);
  (*input->vars)[id + "_z"] = Var(id + "_z", ftot_reduced[2]);
}
