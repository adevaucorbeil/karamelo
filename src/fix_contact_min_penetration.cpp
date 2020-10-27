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

#include "fix_contact_min_penetration.h"
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

FixContactMinPenetration::FixContactMinPenetration(MPM *mpm, vector<string> args)
    : Fix(mpm, args) {
  if (args.size() < 3) {
    error->all(FLERR, "Error: not enough arguments.\n");
  }

  if (args[2].compare("restart") ==
      0) { // If the keyword restart, we are expecting to have read_restart()
           // launched right after.
    solid1 = solid2 = -1;
    mu = 0;
    return;
  }

  if (args.size() < Nargs) {
    error->all(FLERR, "Error: not enough arguments.\n" + usage);
  }

  solid1 = domain->find_solid(args[2]);
  if (solid1 < 0) {
    error->all(FLERR, "Error: solid " + args[2] + " unknown.\n");
  }

  solid2 = domain->find_solid(args[3]);
  if (solid1 < 0) {
    error->all(FLERR, "Error: solid " + args[3] + " unknown.\n");
  }

  cout << "Creating new fix FixContactMinPenetration with ID: " << args[0] << endl;
  id = args[0];
  requires_ghost_particles = true;

  mu = input->parsev(args[4]);
}

FixContactMinPenetration::~FixContactMinPenetration() {}

void FixContactMinPenetration::init() {}

void FixContactMinPenetration::setup() {}

void FixContactMinPenetration::setmask() {
  mask = 0;
  mask |= INITIAL_INTEGRATE;
  //mask |= POST_ADVANCE_PARTICLES;
}
void FixContactMinPenetration::initial_integrate() {
  // cout << "In FixContactMinPenetration::initial_integrate()\n";

  // Go through all the particles in the group and set b to the right value:
  Eigen::Vector3d f, dx;

  Solid *s1, *s2;
  Eigen::Vector3d ftot, ftot_reduced, vtemp1, vtemp2, dv, vt;

  double Rp, Rp1, Rp2, r, inv_r, p, Estar, max_cellsize, vtnorm, fmag;

  ftot.setZero();

  s1 = domain->solids[solid1];
  s2 = domain->solids[solid2];

  Estar = 1.0 / ((1 - s1->mat->nu * s1->mat->nu) / s1->mat->E +
                 (1 - s2->mat->nu * s2->mat->nu) / s2->mat->E);

  max_cellsize = MAX(s1->grid->cellsize, s2->grid->cellsize);

  if (domain->dimension == 2) {
    for (int ip1 = 0; ip1 < s1->np_local; ip1++) {
      for (int ip2 = 0; ip2 < s2->np_local; ip2++) {
        dx = s2->x[ip2] - s1->x[ip1];

        // Extremely gross screening:
        if ((dx[0] < max_cellsize) && (dx[1] < max_cellsize) &&
	    (dx[0] > -max_cellsize) && (dx[1] > -max_cellsize) ) {
	  if (domain->axisymmetric) {
	    Rp1 = 0.5 * sqrt(s1->vol[ip1]/s1->x[ip1][0]);
	    Rp2 = 0.5 * sqrt(s2->vol[ip2]/s2->x[ip2][0]);
	  } else {
	    Rp1 = 0.5 * sqrt(s1->vol[ip1]);
	    Rp2 = 0.5 * sqrt(s2->vol[ip2]);
	  }
          Rp = Rp1 + Rp2;


          // Gross screening:
          if ((dx[0] < Rp) && (dx[1] < Rp) && (dx[0] > -Rp)
	      && (dx[1] > -Rp)) {

            r = dx.norm();

	    // Finer screening:
            if (r < Rp) {
	      inv_r = 1.0/r;
              fmag =
                  s1->mass[ip1] * s2->mass[ip2] /
                  ((s1->mass[ip1] + s2->mass[ip2]) * update->dt * update->dt) *
                  (1 - Rp * inv_r);
              f = fmag * dx;

	      if (mu != 0) {
                dv = s2->v[ip2] - s1->v[ip1];
                vt = dv - dv.dot(dx) * inv_r * inv_r* dx;
		vtnorm = vt.norm();
                if (vtnorm != 0) {
		  vt /= vtnorm;
		  f += mu * fmag * r * vt;
		}
	      }

	      s1->mbp[ip1] += f;
	      s2->mbp[ip2] -= f;
              ftot += f;

              // if ((s1->ptag[ip1] == 12862 && s2->ptag[ip2] == 12902) ||
              //     (s1->ptag[ip1] == 12835 && s2->ptag[ip2] == 12875)) {
              //   cout << "dx[" << s1->ptag[ip1] << "," << s2->ptag[ip2] << "]= ["
              //        << dx[0] << "\t" << dx[1] << "\t" << dx[2] << "]\n";
              //   cout << "f[" << s1->ptag[ip1] << "," << s2->ptag[ip2] << "]= ["
              //        << f[0] << "\t" << f[1] << "\t" << f[2] << "]\n";
              // }
            }
          }
        }
      }
    }
  } else if (domain->dimension == 3) {
    for (int ip1 = 0; ip1 < s1->np_local; ip1++) {
      for (int ip2 = 0; ip2 < s2->np_local; ip2++) {
        dx = s2->x[ip2] - s1->x[ip1];

        // Extremely gross screening:
        if ((dx[0] < max_cellsize) && (dx[1] < max_cellsize) &&
            (dx[2] < max_cellsize) && (dx[0] > -max_cellsize) &&
            (dx[1] > -max_cellsize) && (dx[2] > -max_cellsize)) {
          Rp1 = 0.5 * cbrt(s1->vol[ip1]);
          Rp2 = 0.5 * cbrt(s2->vol[ip2]);
          Rp = Rp1 + Rp2;

          // Gross screening:
          if ((dx[0] < Rp) && (dx[1] < Rp) && (dx[2] < Rp) && (dx[0] > -Rp) &&
              (dx[1] > -Rp) && (dx[2] > -Rp)) {

            r = dx.norm();

	    // Finer screening:
            if (r < Rp) {
	      inv_r = 1.0/r;
              fmag =
                  s1->mass[ip1] * s2->mass[ip2] /
                  ((s1->mass[ip1] + s2->mass[ip2]) * update->dt * update->dt) *
                  (1 - Rp * inv_r);
              f = fmag * dx;

	      if (mu != 0) {
                dv = s2->v[ip2] - s1->v[ip1];
                vt = dv - dv.dot(dx) * inv_r * inv_r* dx;
		vtnorm = vt.norm();
                if (vtnorm != 0) {
		  vt /= vtnorm;
		  f += mu * fmag * r * vt;
		}
	      }
	      s1->mbp[ip1] += f;
	      s2->mbp[ip2] -= f;
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

// void FixContactMinPenetration::post_advance_particles() {
//   // cout << "In FixContactMinPenetration::initial_integrate()\n";

//   // Go through all the particles in the group and set b to the right value:
//   Eigen::Vector3d f, dx;

//   Solid *s1, *s2;
//   Eigen::Vector3d ftot, ftot_reduced, vtemp1, vtemp2;

//   double Rp, Rp1, Rp2, r, p, Estar, max_cellsize;

//   ftot.setZero();

//   s1 = domain->solids[solid1];
//   s2 = domain->solids[solid2];

//   Estar = 1.0 / ((1 - s1->mat->nu * s1->mat->nu) / s1->mat->E +
//                  (1 - s2->mat->nu * s2->mat->nu) / s2->mat->E);

//   max_cellsize = MAX(s1->grid->cellsize, s2->grid->cellsize);

//   if (domain->dimension == 2) {
//     for (int ip1 = 0; ip1 < s1->np_local; ip1++) {
//       for (int ip2 = 0; ip2 < s2->np_local; ip2++) {
//         dx = s2->x[ip2] - s1->x[ip1];

//         // Extremely gross screening:
//         if ((dx[0] < max_cellsize) && (dx[1] < max_cellsize) &&
//             (dx[2] < max_cellsize) && (dx[0] > -max_cellsize) &&
//             (dx[1] > -max_cellsize) && (dx[2] > -max_cellsize)) {
//           Rp1 = 0.5 * sqrt(s1->vol[ip1]);
//           Rp2 = 0.5 * sqrt(s2->vol[ip2]);
//           Rp = Rp1 + Rp2;

//           // Gross screening:
//           if ((dx[0] < Rp) && (dx[1] < Rp) && (dx[2] < Rp) && (dx[0] > -Rp) &&
//               (dx[1] > -Rp) && (dx[2] > -Rp)) {

//             r = dx.norm();
// 	    // Finer screening:
//             if (r < Rp) {
//               f = (1.0 / (s1->mass[ip1] + s2->mass[ip2])) * (1 - Rp / r) * dx;
//               //s1->x[ip1] += s2->mass[ip2] * f;
//               //s2->x[ip2] -= s1->mass[ip1] * f;
// 	      f /= update->dt; 
//               s1->v[ip1] += s2->mass[ip2] * f;
//               s2->v[ip2] -= s1->mass[ip1] * f;
//               ftot += f * (s1->mass[ip1] * s2->mass[ip2])/update->dt;
//             }
//           }
//         }
//       }
//     }
//   } else if (domain->dimension == 3) {
//     for (int ip1 = 0; ip1 < s1->np_local; ip1++) {
//       for (int ip2 = 0; ip2 < s2->np_local; ip2++) {
//         dx = s2->x[ip2] - s1->x[ip1];

//         // Extremely gross screening:
//         if ((dx[0] < max_cellsize) && (dx[1] < max_cellsize) &&
//             (dx[2] < max_cellsize) && (dx[0] > -max_cellsize) &&
//             (dx[1] > -max_cellsize) && (dx[2] > -max_cellsize)) {
//           Rp1 = 0.5 * pow(s1->vol[ip1], 0.333333333);
//           Rp2 = 0.5 * pow(s2->vol[ip2], 0.333333333);
//           Rp = Rp1 + Rp2;

//           // Gross screening:
//           if ((dx[0] < Rp) && (dx[1] < Rp) && (dx[2] < Rp) && (dx[0] > -Rp) &&
//               (dx[1] > -Rp) && (dx[2] > -Rp)) {

//             r = dx.norm();
//             // Finer screening:
//             if (r < Rp) {
//               f = (1.0 / (s1->mass[ip1] + s2->mass[ip2])) * (1 - Rp / r) * dx;
//               //s1->x[ip1] += s2->mass[ip2] * f;
//               //s2->x[ip2] -= s1->mass[ip1] * f;
// 	      f /= update->dt; 
//               s1->v[ip1] += s2->mass[ip2] * f;
//               s2->v[ip2] -= s1->mass[ip1] * f;
//               ftot += f * (s1->mass[ip1] * s2->mass[ip2])/update->dt;
//             }
//           }
//         }
//       }
//     }
//   }

//   // Reduce ftot:
//   MPI_Allreduce(ftot.data(), ftot_reduced.data(), 3, MPI_DOUBLE, MPI_SUM,
//                 universe->uworld);

//   (*input->vars)[id + "_x"] = (*input->vars)[id + "_x"] + Var(id + "_x", ftot_reduced[0]);
//   (*input->vars)[id + "_y"] = (*input->vars)[id + "_y"] + Var(id + "_y", ftot_reduced[1]);
//   (*input->vars)[id + "_z"] = (*input->vars)[id + "_z"] + Var(id + "_z", ftot_reduced[2]);
// }

void FixContactMinPenetration::write_restart(ofstream *of) {
  of->write(reinterpret_cast<const char *>(&solid1), sizeof(int));
  of->write(reinterpret_cast<const char *>(&solid2), sizeof(int));
  of->write(reinterpret_cast<const char *>(&mu), sizeof(double));
}

void FixContactMinPenetration::read_restart(ifstream *ifr) {
  ifr->read(reinterpret_cast<char *>(&solid1), sizeof(int));
  ifr->read(reinterpret_cast<char *>(&solid2), sizeof(int));
  ifr->read(reinterpret_cast<char *>(&mu), sizeof(double));
}
