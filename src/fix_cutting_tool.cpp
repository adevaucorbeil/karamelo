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

#include "fix_cutting_tool.h"
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

FixCuttingTool::FixCuttingTool(MPM *mpm, vector<string> args)
    : Fix(mpm, args) {
  if (args.size() < Nargs) {
    error->all(FLERR, "Error: not enough arguments.\n" + usage);
  }

  int k = 2;
  Kpos = ++k;
  xtpos = ++k;
  ytpos = ++k;
  ztpos = ++k;
  vtxpos = ++k;
  vtypos = ++k;
  vtzpos = ++k;
  xApos = ++k;
  yApos = ++k;
  xBpos = ++k;
  yBpos = ++k;

  if (group->pon[igroup].compare("particles") != 0 &&
      group->pon[igroup].compare("all") != 0) {
    error->all(FLERR, "fix_cuttingtool needs to be given a group of nodes" +
                          group->pon[igroup] + ", " + args[2] +
                          " is a group of " + group->pon[igroup] + ".\n");
  }
  cout << "Creating new fix FixCuttingTool with ID: " << args[0] << endl;
  id = args[0];
}

FixCuttingTool::~FixCuttingTool() {}

void FixCuttingTool::init() {}

void FixCuttingTool::setup() {}

void FixCuttingTool::setmask() {
  mask = 0;
  mask |= INITIAL_INTEGRATE;
}

void FixCuttingTool::initial_integrate() {
  // cout << "In FixCuttingTool::initial_integrate()\n";

  // Go through all the particles in the group and set b to the right value:
  Eigen::Vector3d f;

  int solid = group->solid[igroup];

  Solid *s;
  Eigen::Vector3d ftot, ftot_reduced, n1, n2, n;

  double K = input->parsev(args[Kpos]).result(mpm);

  Eigen::Vector3d xt(input->parsev(args[xtpos]).result(mpm),
                     input->parsev(args[ytpos]).result(mpm),
                     input->parsev(args[ztpos]).result(mpm));
  Eigen::Vector3d vt(input->parsev(args[vtxpos]).result(mpm),
                     input->parsev(args[vtypos]).result(mpm),
                     input->parsev(args[vtzpos]).result(mpm));
  Eigen::Vector3d xA(input->parsev(args[xApos]).result(mpm),
                     input->parsev(args[yApos]).result(mpm),
                     input->parsev(args[ztpos]).result(mpm));
  Eigen::Vector3d xB(input->parsev(args[xBpos]).result(mpm),
                     input->parsev(args[yBpos]).result(mpm),
                     input->parsev(args[ztpos]).result(mpm));

  // The equation of line 1 is: (yA - yt) * x - (xA - xt) * y + yt * xA - yA *
  // xt = 0 The equation of line 2 is: (yB - yt) * x - (xB - xt) * y + yt * xB -
  // yB * xt = 0

  double line1[4], line2[4];
  line1[0] = xA[1] - xt[1];
  line1[1] = xA[0] - xt[0];
  line1[2] = xt[1] * xA[0] - xt[0] * xA[1];
  line1[3] = 1.0 / sqrt(line1[0] * line1[0] + line1[1] * line1[1]);

  if (line1[0] * xB[0] + line1[1] * xB[1] + line1[2] < 0) {
    line1[0] *= -1;
    line1[1] *= -1;
    line1[2] *= -1;
  }

  if (line2[0] * xA[0] + line2[1] * xA[1] + line2[2] < 0) {
    line2[0] *= -1;
    line2[1] *= -1;
    line2[2] *= -1;
  }

  n1[0] = line1[0];
  n1[1] = line1[1];
  n1[2] = 0;
  n1 *= -line1[3];

  n2[0] = line2[0];
  n2[1] = line2[1];
  n2[2] = 0;
  n2 *= -line2[3];

  ftot.setZero();

  double r, p, p1, p2, c1p, c2p, fmag;
  if (solid == -1) {
    for (int isolid = 0; isolid < domain->solids.size(); isolid++) {
      s = domain->solids[isolid];

      if (domain->dimension == 2) {
        for (int ip = 0; ip < s->np_local; ip++) {
          if (s->mass[ip] > 0) {
            if (s->mask[ip] & groupbit) {

              c1p = line1[0] * s->x[ip][0] + line1[1] * s->x[ip][1] + line1[1];
              c2p = line2[0] * s->x[ip][0] + line2[1] * s->x[ip][1] + line2[1];

              if (c1p >= 0 && c2p >= 0) {
                // The particle is inside the tool
                p1 = fabs(c1p * line1[3]);
                p2 = fabs(c2p * line2[3]);

                if (p1 < p2) {
                  p = p1;
                  n = n1;
                } else {
                  p = p2;
                  n = n2;
                }
                fmag = K * s->mat->G * p;

                f = fmag * n;
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

      if (domain->dimension == 3) {
        // Not supported
        error->one(FLERR, "fix_cuttingtool not supported in 3D\n");
      }
    }
  } else {
    s = domain->solids[solid];

    if (domain->dimension == 2) {
      for (int ip = 0; ip < s->np_local; ip++) {
        if (s->mass[ip] > 0) {
          if (s->mask[ip] & groupbit) {

            c1p = line1[0] * s->x[ip][0] + line1[1] * s->x[ip][1] + line1[1];
            c2p = line2[0] * s->x[ip][0] + line2[1] * s->x[ip][1] + line2[1];

            if (c1p >= 0 && c2p >= 0) {
              // The particle is inside the tool
              p1 = fabs(c1p * line1[3]);
              p2 = fabs(c2p * line2[3]);

              if (p1 < p2) {
                p = p1;
                n = n1;
              } else {
                p = p2;
                n = n2;
              }
              fmag = K * s->mat->G * p;

              f = fmag * n;
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

    if (domain->dimension == 3) {
      // Not supported
      error->one(FLERR, "fix_cuttingtool not supported in 3D\n");
    }
  }
  // Reduce ftot:
  MPI_Allreduce(ftot.data(), ftot_reduced.data(), 3, MPI_DOUBLE, MPI_SUM,
		universe->uworld);

  (*input->vars)[id + "_x"] = Var(id + "_x", ftot_reduced[0]);
  (*input->vars)[id + "_y"] = Var(id + "_y", ftot_reduced[1]);
  (*input->vars)[id + "_z"] = Var(id + "_z", ftot_reduced[2]);
}
