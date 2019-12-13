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

#include "fix_indent_pinball.h"
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

#define HALF_SQRT_2 0.7071067811865476

FixIndentPinball::FixIndentPinball(MPM *mpm, vector<string> args)
    : Fix(mpm, args) {
  if (args.size() < Nargs) {
    error->all(FLERR, "Error: not enough arguments.\n" + usage);
  }

  type_pos = 3;
  Kpos = 4;
  Rpos = 5;
  xpos = 6;
  ypos = 7;
  zpos = 8;
  vxpos = 9;
  vypos = 10;
  vzpos = 11;

  if (group->pon[igroup].compare("particles") != 0 &&
      group->pon[igroup].compare("all") != 0) {
    error->all(FLERR, "fix_indent_pinball needs to be given a group of nodes" +
                          group->pon[igroup] + ", " + args[2] +
                          " is a group of " + group->pon[igroup] + ".\n");
  }
  cout << "Creating new fix FixIndentPinball with ID: " << args[0] << endl;
  id = args[0];

  type = args[type_pos];
  if (args[type_pos].compare("sphere") == 0) {
    type = "sphere";
  } else {
    error->all(FLERR, "Error indent type " + args[type_pos] +
                          " unknown. Only type sphere is supported.\n");
  }
}

FixIndentPinball::~FixIndentPinball() {}

void FixIndentPinball::init() {}

void FixIndentPinball::setup() {}

void FixIndentPinball::setmask() {
  mask = 0;
  mask |= INITIAL_INTEGRATE;
}

void FixIndentPinball::initial_integrate() {
  // cout << "In FixIndentPinball::initial_integrate()\n";

  // Go through all the particles in the group and set b to the right value:
  Eigen::Vector3d f;

  int solid = group->solid[igroup];

  Solid *s;
  Eigen::Vector3d ftot, ftot_reduced;

  double K = input->parsev(args[Kpos]).result(mpm);
  double R = input->parsev(args[Rpos]).result(mpm);
  Eigen::Vector3d xs(input->parsev(args[xpos]).result(mpm),
                     input->parsev(args[ypos]).result(mpm),
                     input->parsev(args[zpos]).result(mpm));
  Eigen::Vector3d vs(input->parsev(args[vxpos]).result(mpm),
                     input->parsev(args[vypos]).result(mpm),
                     input->parsev(args[vzpos]).result(mpm));
  Eigen::Vector3d xsp;

  double Rs, Rp, r, p, pdot, f1, f2, fmag;

  ftot.setZero();

  if (solid == -1) {
    for (int isolid = 0; isolid < domain->solids.size(); isolid++) {
      s = domain->solids[isolid];

      for (int ip = 0; ip < s->np_local; ip++) {
        if (s->mass[ip] > 0) {
          if (s->mask[ip] & groupbit) {

            // Gross screening:
            xsp = s->x[ip] - xs;

            Rp = HALF_SQRT_2 * s->grid->cellsize;
            Rs = R + Rp;

            if ((xsp[0] < Rs) && (xsp[1] < Rs) && (xsp[2] < Rs) &&
                (xsp[0] > -Rs) && (xsp[1] > -Rs) && (xsp[2] > -Rs)) {

              r = xsp.norm();
              // Finer screening:
              if (r < Rs) {
                p = Rs - r; // penetration

                pdot = xsp.dot(vs - s->v[ip]); // Not yet normalized

                if (pdot > 0) {
                  pdot /= r; // Normalized only here to save computation time
                  f1 = s->rho[ip] * Rp * Rp * Rp * pdot / update->dt;
                } else {
                  f1 = 0;
                }

                f2 = s->mat->G * sqrt(Rp * p * p * p);

                fmag = K * MIN(f1, f2);

                f = fmag * xsp / r;
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

    for (int ip = 0; ip < s->np_local; ip++) {
      if (s->mass[ip] > 0) {
        if (s->mask[ip] & groupbit) {

          // Gross screening:
          xsp = s->x[ip] - xs;

          Rp = HALF_SQRT_2 * s->grid->cellsize;
          Rs = R + Rp;

          if ((xsp[0] < Rs) && (xsp[1] < Rs) && (xsp[2] < Rs) &&
              (xsp[0] > -Rs) && (xsp[1] > -Rs) && (xsp[2] > -Rs)) {

            r = xsp.norm();
            // Finer screening:
            if (r < Rs) {
              p = Rs - r; // penetration

              pdot = xsp.dot(vs - s->v[ip]); // Not yet normalized

              if (pdot > 0) {
                pdot /= r; // Normalized only here to save computation time
                f1 = s->rho[ip] * Rp * Rp * Rp * pdot / update->dt;
              } else {
                f1 = 0;
              }

              f2 = s->mat->G * sqrt(Rp * p * p * p);

              fmag = K * MIN(f1, f2);

              f = fmag * xsp / r;
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
