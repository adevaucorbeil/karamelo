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

#include "fix_indent_hertz_nodes_symz.h"
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

FixIndentHertzNodesSymZ::FixIndentHertzNodesSymZ(MPM *mpm, vector<string> args)
    : Fix(mpm, args) {
  if (args.size() < Nargs) {
    error->all(FLERR, "Error: not enough arguments.\n" + usage);
  }

  type_pos = 3;
  Rpos = 4;
  xpos = 5;
  ypos = 6;
  zpos = 7;
  vxpos = 8;
  vypos = 9;
  vzpos = 10;

  if (group->pon[igroup].compare("nodes") != 0 &&
      group->pon[igroup].compare("all") != 0) {
    error->all(FLERR, "fix_indent_hertz needs to be given a group of nodes" +
                          group->pon[igroup] + ", " + args[2] +
                          " is a group of " + group->pon[igroup] + ".\n");
  }
  cout << "Creating new fix FixIndentHertzNodesSymZ with ID: " << args[0] << endl;
  id = args[0];

  type = args[type_pos];
  if (args[type_pos].compare("sphere") == 0) {
    type = "sphere";
  } else {
    error->all(FLERR, "Error indent type " + args[type_pos] +
                          " unknown. Only type sphere is supported.\n");
  }
}

FixIndentHertzNodesSymZ::~FixIndentHertzNodesSymZ() {}

void FixIndentHertzNodesSymZ::init() {}

void FixIndentHertzNodesSymZ::setup() {}

void FixIndentHertzNodesSymZ::setmask() {
  mask = 0;
  mask |= POST_PARTICLES_TO_GRID;
}

void FixIndentHertzNodesSymZ::post_particles_to_grid() {
  // cout << "In FixIndentHertzNodesSymZ::post_particles_to_grid()\n";

  // Go through all the particles in the group and set b to the right value:
  Eigen::Vector3d f;

  int solid = group->solid[igroup];

  Grid *g;
  Solid *s;
  Eigen::Vector3d ftot, ftot_reduced;

  double Rs = input->parsev(args[Rpos]).result(mpm);
  Eigen::Vector3d xs(input->parsev(args[xpos]).result(mpm),
                     input->parsev(args[ypos]).result(mpm),
                     input->parsev(args[zpos]).result(mpm));
  Eigen::Vector3d vs(input->parsev(args[vxpos]).result(mpm),
                     input->parsev(args[vypos]).result(mpm),
                     input->parsev(args[vzpos]).result(mpm));
  Eigen::Vector3d xsn;

  double r, p, pdot, f1, f2, fmag, k;

  int n, n_reduced;
  ftot.setZero();
  n = 0;

  if (solid == -1) {
    for (int isolid = 0; isolid < domain->solids.size(); isolid++) {
      s = domain->solids[isolid];
      g = domain->solids[isolid]->grid;

      if (domain->dimension == 2) {

	k = 0.25 * M_PI * s->mat->E / (1 - s->mat->nu * s->mat->nu);

        for (int in = 0; in < g->nnodes_local + g->nnodes_ghost; in++) {
          if (g->mass[in] > 0) {
            if (g->mask[in] & groupbit) {

              // Gross screening:
              xsn = g->x[in] - xs;

              if ((xsn[0] < Rs) && (xsn[1] < Rs) && (xsn[0] > -Rs) &&
                  (xsn[1] > -Rs)) {

                r = xsn.norm();

                // Finer screening:
                if (r < Rs) {
                  p = Rs - r; // penetration

                  if (p > 0) {
		    n++;
                    fmag = k * p;

                    // fmag = K * MIN(f1, f2);
		    if (g->x[in][1] < 1e-12 && g->x[in][1] > -1e-12)
		      fmag /= 2;
                    f = fmag * xsn / r;
                    g->mb[in] += f;
                    if (in < g->nnodes_local) ftot += f;
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

	k = four_thirds * s->mat->E / (1 - s->mat->nu * s->mat->nu);

        for (int in = 0; in < g->nnodes_local + g->nnodes_ghost; in++) {
          if (g->mass[in] > 0) {
            if (g->mask[in] & groupbit) {

              // Gross screening:
              xsn = g->x[in] - xs;

              if ((xsn[0] < Rs) && (xsn[1] < Rs) && (xsn[2] < Rs) &&
                  (xsn[0] > -Rs) && (xsn[1] > -Rs) && (xsn[2] > -Rs)) {

                r = xsn.norm();
                // Finer screening:
                if (r < Rs) {
                  p = Rs - r; // penetration

                  if (p > 1.0e-10) {
		    n++;
                    fmag = k * sqrt(Rs * p * p * p);

                    // fmag = K * MIN(f1, f2);

		    if (g->x[in][1] < 1e-12 && g->x[in][1] > -1e-12)
		      fmag /= 2;
                    f = fmag * xsn / r;
                    g->mb[in] += f;
                    if (in < g->nnodes_local) ftot += f;
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
    g = domain->solids[solid]->grid;

    if (domain->dimension == 2) {

      k = 0.25 * M_PI * s->mat->E / (1 - s->mat->nu * s->mat->nu);

      for (int in = 0; in < g->nnodes_local + g->nnodes_ghost; in++) {
        if (g->mass[in] > 0) {
          if (g->mask[in] & groupbit) {

            // Gross screening:
            xsn = g->x[in] - xs;

            if ((xsn[0] < Rs) && (xsn[1] < Rs) && (xsn[0] > -Rs) &&
                (xsn[1] > -Rs)) {

              r = xsn.norm();
              // Finer screening:
              if (r < Rs) {
                p = Rs - r; // penetration

                if (p > 0) {
		  n++;
		  fmag = k * p;

                  // fmag = K * MIN(f1, f2);

		  if (g->x[in][1] < 1e-12 && g->x[in][1] > -1e-12)
		    fmag /= 2;
                  f = fmag * xsn / r;
                  g->mb[in] += f;
                  if (in < g->nnodes_local) ftot += f;
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

      k = four_thirds * s->mat->E / (1 - s->mat->nu * s->mat->nu);
      for (int in = 0; in < g->nnodes_local + g->nnodes_ghost; in++) {
        if (g->mass[in] > 0) {
          if (g->mask[in] & groupbit) {

            // Gross screening:
            xsn = g->x[in] - xs;

            if ((xsn[0] < Rs) && (xsn[1] < Rs) && (xsn[2] < Rs) &&
                (xsn[0] > -Rs) && (xsn[1] > -Rs) && (xsn[2] > -Rs)) {

              r = xsn.norm();
              // Finer screening:
              if (r < Rs) {
                p = Rs - r; // penetration

                if (p > 0) {
		  n++;
		  fmag = k * sqrt(Rs * p * p * p);
                  // fmag = K * MIN(f1, f2);

		  if (g->x[in][1] < 1e-12 && g->x[in][1] > -1e-12)
		    fmag /= 2 ;
                  f = fmag * xsn / r;
                  g->mb[in] += f;
                  if (in < g->nnodes_local) ftot += f;
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
