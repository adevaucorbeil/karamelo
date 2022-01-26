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

#include <fix_indent_minimize_penetration.h>
#include <domain.h>
#include <error.h>
#include <group.h>
#include <input.h>
#include <solid.h>
#include <universe.h>
#include <update.h>
#include <matrix.h>
#include <iostream>
#include <string>
#include <vector>

using namespace std;
using namespace FixConst;


#define four_thirds 1.333333333

FixIndentMinimizePenetration::FixIndentMinimizePenetration(MPM *mpm, vector<string> args)
    : Fix(mpm, args, INITIAL_INTEGRATE) {
  if (args.size() < 3) {
    error->all(FLERR, "Error: not enough arguments.\n");
  }

  if (args[2].compare("restart") ==
      0) { // If the keyword restart, we are expecting to have read_restart()
           // launched right after.
    igroup = stoi(args[3]);
    if (igroup == -1 && universe->me == 0) {
      cout << "Could not find group number " << args[3] << endl;
    }
    groupbit = group->bitmask[igroup];

    R = mu = 0;
    return;
  }

  if (args.size() < Nargs) {
    error->all(FLERR, "Error: not enough arguments.\n" + usage);
  }
  if (args.size() > Nargs) {
    error->all(FLERR, "Error: too many arguments.\n" + usage);
  }

  if (group->pon[igroup].compare("particles") != 0 &&
      group->pon[igroup].compare("all") != 0) {
    error->all(FLERR, "fix_indent_hertz needs to be given a group of particles" +
                          group->pon[igroup] + ", " + args[2] +
                          " is a group of " + group->pon[igroup] + ".\n");
  }
  if (universe->me == 0) {
    cout << "Creating new fix FixIndentMinimizePenetration with ID: " << args[0] << endl;
  }
  id = args[0];

  type = args[3];
  if (args[3].compare("sphere") == 0) {
    type = "sphere";
  } else {
    error->all(FLERR, "Error indent type " + args[3] +
                          " unknown. Only type sphere is supported.\n");
  }

  R = input->parsev(args[4]);
  xvalue = input->parsev(args[5]);
  yvalue = input->parsev(args[6]);
  zvalue = input->parsev(args[7]);
  vxvalue = input->parsev(args[8]);
  vyvalue = input->parsev(args[9]);
  vzvalue = input->parsev(args[10]);
  mu = input->parsev(args[11]);
}

void FixIndentMinimizePenetration::initial_integrate() {
  // Go through all the particles in the group and set b to the right value:
  Vector3d f;

  int solid = group->solid[igroup];

  Solid *s;
  Vector3d ftot, ftot_reduced, ffric;

  Vector3d xs(xvalue.result(mpm),
                     yvalue.result(mpm),
                     zvalue.result(mpm));
  Vector3d vs(vxvalue.result(mpm),
                     vyvalue.result(mpm),
                     vzvalue.result(mpm));
  Vector3d xsp, vps, vt;
  Vector3d ey = {0, 1, 0};

  double Rs, Rp, r, p, fmag, vtnorm, vndotxsp;

  double A, A_reduced, cellsizeSq;
  ftot = Vector3d();
  A = 0;

  if (solid == -1) {
    for (int isolid = 0; isolid < domain->solids.size(); isolid++) {
      s = domain->solids[isolid];
      cellsizeSq = s->grid->cellsize * s->grid->cellsize;

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

		    A += cellsizeSq * s->F[ip](1,1) * s->F[ip](2,2);
		    s->mbp[ip] += f;
                    ftot += f;
                  } else {
                    fmag = 0;
                    f = Vector3d();
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
			// cout << "Particle " << s->ptag[ip]
			//      << " f_friction=[" << ffric[0] <<"," << ffric[1] <<"," << ffric[2] << "] "
			//      << "f=[" << f[0] <<"," << f[1] <<"," << f[2] << "] "
			//      << "vps=[" << vps[0] <<"," << vps[1] <<"," << vps[2] << "] "
			//       << "vt=[" << vt[0] <<"," << vt[1] <<"," << vt[2] << "] "
			// //      << "vs=[" << vs[0] <<"," << vs[1] <<"," << vs[2] << "] "
			// //      << "vp=[" << s->v[ip][0] <<"," << s->v[ip][1] <<"," << s->v[ip][2] << "] "
			// // //      << "ap=[" << s->a[ip][0] <<"," << s->a[ip][1] <<"," << s->a[ip][2] << "] "
			//      << "xsp=[" << xsp[0] <<"," << xsp[1] <<"," << xsp[2] << "] "
			//      << "dt=" << update->dt << endl;
			// // " mass=" << s->mass[ip] << " vtnorm=" << vtnorm << "\n";
			// // //   cout << endl;

			f += ffric;
		      }
		    }

		    A += cellsizeSq * s->F[ip](0,0) * s->F[ip](2,2);
		    // cout << "Particle " << s->ptag[ip] << " h^2=" << cellsizeSq << " F11=" << s->F[ip](0,0) << " F33=" << s->F[ip](2,2) << " dA=" << cellsizeSq * s->F[ip](0,0) * s->F[ip](2,2) << endl;
		    s->mbp[ip] += f;
                    ftot += f;
                  } else {
                    fmag = 0;
                    f = Vector3d();
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
    cellsizeSq = s->grid->cellsize * s->grid->cellsize;

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

		  A += cellsizeSq * s->F[ip](0,0) * s->F[ip](2,2);
		  s->mbp[ip] += f;
		  ftot += f;
		} else {
		  fmag = 0;
		  f = Vector3d();
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

		  A += cellsizeSq * s->F[ip](0,0) * s->F[ip](2,2);
		  s->mbp[ip] += f;
		  ftot += f;
		} else {
		  fmag = 0;
		  f = Vector3d();
		}
	      }
	    }
	  }
	}
      }
    }
  }

  // Reduce ftot:
  MPI_Allreduce(&A, &A_reduced, 1, MPI_DOUBLE, MPI_SUM, universe->uworld);
  MPI_Allreduce(ftot.elements, ftot_reduced.elements, 3, MPI_DOUBLE, MPI_SUM,
                universe->uworld);

  (*input->vars)[id + "_s"] = Var(id + "_s", A_reduced);
  (*input->vars)[id + "_x"] = Var(id + "_x", ftot_reduced[0]);
  (*input->vars)[id + "_y"] = Var(id + "_y", ftot_reduced[1]);
  (*input->vars)[id + "_z"] = Var(id + "_z", ftot_reduced[2]);
}

void FixIndentMinimizePenetration::write_restart(ofstream *of) {
  of->write(reinterpret_cast<const char *>(&R), sizeof(double));
  of->write(reinterpret_cast<const char *>(&mu), sizeof(double));
  xvalue.write_to_restart(of);
  yvalue.write_to_restart(of);
  zvalue.write_to_restart(of);

  vxvalue.write_to_restart(of);
  vyvalue.write_to_restart(of);
  vzvalue.write_to_restart(of);
}

void FixIndentMinimizePenetration::read_restart(ifstream *ifr) {
  ifr->read(reinterpret_cast<char *>(&R), sizeof(double));
  ifr->read(reinterpret_cast<char *>(&mu), sizeof(double));
  xvalue.read_from_restart(ifr);
  yvalue.read_from_restart(ifr);
  zvalue.read_from_restart(ifr);

  vxvalue.read_from_restart(ifr);
  vyvalue.read_from_restart(ifr);
  vzvalue.read_from_restart(ifr);
  type = "sphere";
}
