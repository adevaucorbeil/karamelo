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

#include "fix_impenetrable_surface.h"
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

FixImpenetrableSurface::FixImpenetrableSurface(MPM *mpm, vector<string> args)
    : Fix(mpm, args) {
  if (args.size() < Nargs) {
    error->all(FLERR, "Error: not enough arguments.\n" + usage);
  }

  int k = 2;
  Kpos = ++k;
  xspos = ++k;
  yspos = ++k;
  zspos = ++k;
  nxpos = ++k;
  nypos = ++k;
  nzpos = ++k;

  if (group->pon[igroup].compare("particles") != 0 &&
      group->pon[igroup].compare("all") != 0) {
    error->all(FLERR,
               "fix_impenetrablesurface needs to be given a group of nodes" +
                   group->pon[igroup] + ", " + args[2] + " is a group of " +
                   group->pon[igroup] + ".\n");
  }
  cout << "Creating new fix FixImpenetrableSurface with ID: " << args[0]
       << endl;
  id = args[0];
}

FixImpenetrableSurface::~FixImpenetrableSurface() {}

void FixImpenetrableSurface::init() {}

void FixImpenetrableSurface::setup() {}

void FixImpenetrableSurface::setmask() {
  mask = 0;
  mask |= INITIAL_INTEGRATE;
}

void FixImpenetrableSurface::initial_integrate() {
  // cout << "In FixImpenetrableSurface::initial_integrate()\n";

  // Go through all the particles in the group and set b to the right value:
  Eigen::Vector3d f;

  int solid = group->solid[igroup];

  Solid *s;
  Eigen::Vector3d ftot, ftot_reduced;

  double K = input->parsev(args[Kpos]).result(mpm);

  Eigen::Vector3d xs(input->parsev(args[xspos]).result(mpm),
                     input->parsev(args[yspos]).result(mpm),
                     input->parsev(args[zspos]).result(mpm));
  Eigen::Vector3d n(input->parsev(args[nxpos]).result(mpm),
                    input->parsev(args[nypos]).result(mpm),
                    input->parsev(args[nzpos]).result(mpm)); // Outgoing normal

  // Normalize n:
  n /= n.norm();
  ftot.setZero();

  double p, fmag;

  double D = -n[0] * xs[0] - n[1] * xs[1] - n[2] * xs[2];

  // cout << "line 1: " << line1[0] << "x + " << line1[1] << "y + " << line1[2]
  // << endl; cout << "line 2: " << line2[0] << "x + " << line2[1] << "y + " <<
  // line2[2] << endl;

  if (solid == -1) {
    for (int isolid = 0; isolid < domain->solids.size(); isolid++) {
      s = domain->solids[isolid];

      for (int ip = 0; ip < s->np_local; ip++) {
        if (s->mass[ip] > 0) {
          if (s->mask[ip] & groupbit) {

            p = -(n[0] * s->x[ip][0] + n[1] * s->x[ip][1] + n[2] * s->x[ip][2] +
                  D);
            // if (s->ptag[ip] == 1) {
            //   cout << id << "- Particle " << s->ptag[ip] << "\t";
            //   cout << "p = " << p << "\t";
            //   cout << "xp = [" << s->x[ip][0] << "," << s->x[ip][1] << ","
            //        << s->x[ip][2] << "]\t" << endl;
            //   cout << "xs = [" << xs[0] << "," << xs[1] << "," << xs[2] << "]\t"
            //        << endl;
            //   cout << "n = [" << n[0] << "," << n[1] << "," << n[2] << "]"
            //        << endl;
            // }

            if (p >= 0) {
              // cout << "Particle " << s->ptag[ip] << " is inside " << id << "\n";
              fmag = K * s->mat->G * p * (1.0 - s->damage[ip]);

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
  } else {
    s = domain->solids[solid];

    for (int ip = 0; ip < s->np_local; ip++) {
      if (s->mass[ip] > 0) {
        if (s->mask[ip] & groupbit) {

          p = -(n[0] * s->x[ip][0] + n[1] * s->x[ip][1] + n[2] * s->x[ip][2] +
                D);
          // if (s->ptag[ip] == 1) {
          //   cout << id << "- Particle " << s->ptag[ip] << "\t";
          //   cout << "p = " << p << "\t";
          //   cout << "xp = [" << s->x[ip][0] << "," << s->x[ip][1] << ","
          //        << s->x[ip][2] << "]\t" << endl;
          //   cout << "xs = [" << xs[0] << "," << xs[1] << "," << xs[2] << "]\t"
          //        << endl;
          //   cout << "n = [" << n[0] << "," << n[1] << "," << n[2] << "]"
          //        << endl;
          // }

          if (p >= 0) {
            // cout << "Particle " << s->ptag[ip] << " is inside\n";
            fmag = K * s->mat->G * p * (1.0 - s->damage[ip]);

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
  // Reduce ftot:
  MPI_Allreduce(ftot.data(), ftot_reduced.data(), 3, MPI_DOUBLE, MPI_SUM,
                universe->uworld);

  (*input->vars)[id + "_x"] = Var(id + "_x", ftot_reduced[0]);
  (*input->vars)[id + "_y"] = Var(id + "_y", ftot_reduced[1]);
  (*input->vars)[id + "_z"] = Var(id + "_z", ftot_reduced[2]);
}
