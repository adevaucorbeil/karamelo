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

#include "fix_contact_hertz.h"
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

FixContactHertz::FixContactHertz(MPM *mpm, vector<string> args)
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

  cout << "Creating new fix FixContactHertz with ID: " << args[0] << endl;
  id = args[0];
  requires_ghost_particles = true;
}

FixContactHertz::~FixContactHertz() {}

void FixContactHertz::init() {}

void FixContactHertz::setup() {}

void FixContactHertz::setmask() {
  mask = 0;
  mask |= INITIAL_INTEGRATE;
}

void FixContactHertz::initial_integrate() {
  // cout << "In FixContactHertz::initial_integrate()\n";

  // Go through all the particles in the group and set b to the right value:
  Eigen::Vector3d f, dx;

  Solid *s1, *s2;
  Eigen::Vector3d ftot, ftot_reduced;

  double Rp, Rp1, Rp2, r, p, fmag, Estar;

  ftot.setZero();

  s1 = domain->solids[solid1];
  s2 = domain->solids[solid2];

  Estar = 1.0 / ((1 - s1->mat->nu * s1->mat->nu) / s1->mat->E +
                 (1 - s2->mat->nu * s2->mat->nu) / s2->mat->E);

  if (domain->dimension == 3) {
    for (int ip1 = 0; ip1 < s1->np_local; ip1++) {
      for (int ip2 = 0; ip2 < s2->np_local; ip2++) {
        dx = s2->x[ip2] - s1->x[ip1];
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

            fmag =
                four_thirds * Estar * sqrt(Rp1 * Rp2 / (Rp1 + Rp2) * p * p * p);

            f = fmag * dx / r;
            s1->mbp[ip1] -= f;
            s2->mbp[ip2] += f;
            ftot += f;
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
