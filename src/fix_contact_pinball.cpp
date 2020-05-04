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

#include "fix_contact_pinball.h"
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

FixContactPinball::FixContactPinball(MPM *mpm, vector<string> args)
    : Fix(mpm, args) {
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

  cout << "Creating new fix FixContactPinball with ID: " << args[0] << endl;
  id = args[0];
  requires_ghost_particles = true;

  K = input->parsev(args[4]);
}

FixContactPinball::~FixContactPinball() {}

void FixContactPinball::init() {}

void FixContactPinball::setup() {}

void FixContactPinball::setmask() {
  mask = 0;
  mask |= INITIAL_INTEGRATE;
}

void FixContactPinball::initial_integrate() {
  // cout << "In FixContactPinball::initial_integrate()\n";

  // Go through all the particles in the group and set b to the right value:
  Eigen::Vector3d f, dx, dv;

  Solid *s1, *s2;
  Eigen::Vector3d ftot, ftot_reduced;

  double Rp, Rp1, Rp2, r, p, pdot, fmag1, fmag2, Gstar, max_cellsize;

  ftot.setZero();

  s1 = domain->solids[solid1];
  s2 = domain->solids[solid2];

  Gstar = 1.0 / ((1 - s1->mat->nu * s1->mat->nu) / s1->mat->G +
                 (1 - s2->mat->nu * s2->mat->nu) / s2->mat->G);

  max_cellsize = MAX(s1->grid->cellsize, s2->grid->cellsize);

  for (int ip1 = 0; ip1 < s1->np_local; ip1++) {
    for (int ip2 = 0; ip2 < s2->np_local; ip2++) {
      dx = s2->x[ip2] - s1->x[ip1];

      // Extremely gross screening:
      if ((dx[0] < max_cellsize) && (dx[1] < max_cellsize) &&
          (dx[2] < max_cellsize) && (dx[0] > -max_cellsize) &&
          (dx[1] > -max_cellsize) && (dx[2] > -max_cellsize)) {
        Rp1 = 0.5 * pow(s1->vol[ip1], 0.333333333);
        Rp2 = 0.5 * pow(s2->vol[ip2], 0.333333333);
        Rp = Rp1 + Rp2;

        // Gross screening:
        if ((dx[0] < Rp) && (dx[1] < Rp) && (dx[2] < Rp) && (dx[0] > -Rp) &&
            (dx[1] > -Rp) && (dx[2] > -Rp)) {

          r = dx.norm();

          // Finer screening:
          if (r < Rp) {

            p = Rp - r; // penetration

            dv = s2->v[ip2] - s1->v[ip1];
            pdot = -dx.dot(dv) / r; // penetration speed

            if (pdot > 0) {
              fmag1 = s1->rho[ip1] * s2->rho[ip2] * Rp1 * Rp1 * Rp1 * Rp2 *
                      Rp2 * Rp2 * pdot /
                      ((s1->rho[ip1] * Rp1 * Rp1 * Rp1 +
                        s2->rho[ip2] * Rp2 * Rp2 * Rp2) *
                       update->dt);
              fmag2 = Gstar * sqrt(Rp1 * Rp2 / (Rp1 + Rp2) * p * p * p);

              f = MIN(fmag1, fmag2) * dx / r;
              ftot += f;
	      s1->mbp[ip1] -= f;
              s2->mbp[ip2] += f;
              //s1->f[ip1] -= f;
              //s2->f[ip2] += f;
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
