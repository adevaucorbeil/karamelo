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

#include <fix_contact_hertz.h>
#include <domain.h>
#include <error.h>
#include <group.h>
#include <input.h>
#include <solid.h>
#include <universe.h>
#include <update.h>


using namespace std;
using namespace FixConst;


#define four_thirds 1.333333333

FixContactHertz::FixContactHertz(MPM *mpm, vector<string> args)
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

    solid1 = solid2 = -1;
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

  if (universe->me == 0) {
    cout << "Creating new fix FixContactHertz with ID: " << args[0] << endl;
  }
  id = args[0];
  requires_ghost_particles = true;
}

void FixContactHertz::prepare()
{
  ftot = Vector3d();
}

void FixContactHertz::reduce()
{
  Vector3d ftot_reduced;

  // Reduce ftot:
  MPI_Allreduce(ftot.elements, ftot_reduced.elements, 3, MPI_DOUBLE, MPI_SUM,
                universe->uworld);

  (*input->vars)[id + "_x"] = Var(id + "_x", ftot_reduced[0]);
  (*input->vars)[id + "_y"] = Var(id + "_y", ftot_reduced[1]);
  (*input->vars)[id + "_z"] = Var(id + "_z", ftot_reduced[2]);
}

void FixContactHertz::initial_integrate() {
  // cout << "In FixContactHertz::initial_integrate()\n";

  // Go through all the particles in the group and set b to the right value:
  Vector3d f, dx;

  Solid *s1, *s2;
  Vector3d vtemp1, vtemp2;

  double Rp, Rp1, Rp2, r, p, fmag, Estar, max_cellsize;

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
        (dx[2] < max_cellsize) && (dx[0] > -max_cellsize) &&
            (dx[1] > -max_cellsize) && (dx[2] > -max_cellsize)) {
      Rp1 = 0.5 * sqrt(s1->vol[ip1]);
      Rp2 = 0.5 * sqrt(s2->vol[ip2]);
      Rp = Rp1 + Rp2;

      // Gross screening:
      if ((dx[0] < Rp) && (dx[1] < Rp) && (dx[2] < Rp) && (dx[0] > -Rp) &&
          (dx[1] > -Rp) && (dx[2] > -Rp)) {

        r = dx.norm();

        // Finer screening:
        if (r < Rp) {

              p = Rp - r; // penetration

              fmag = 0.25 * M_PI * Estar *
                     sqrt(Rp1 * Rp2 / (Rp1 + Rp2) * p * p * p);

              f = fmag * dx / r;
              ftot += f;
          s1->mbp[ip1] -= f;
          s2->mbp[ip2] += f;
              // vtemp1 = -update->dt * f / s1->mass[ip1];
              // s1->v[ip1] += vtemp1;
              // s1->x[ip1] += update->dt * vtemp1;

              // vtemp2 = update->dt * f / s2->mass[ip2];
              // s2->v[ip2] += vtemp2;
              // s2->x[ip2] += update->dt * vtemp2;
            }
          }
        }
      }
    }
  }
  else if (domain->dimension == 3) {
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

              fmag = four_thirds * Estar *
                     sqrt(Rp1 * Rp2 / (Rp1 + Rp2) * p * p * p);

              f = fmag * dx / r;
          s1->mbp[ip1] -= f;
          s2->mbp[ip2] += f;
              ftot += f;
              // vtemp1 = -update->dt * f / s1->mass[ip1];
              // s1->v[ip1] += vtemp1;
              // s1->x[ip1] += update->dt * vtemp1;

              // vtemp2 = update->dt * f / s2->mass[ip2];
              // s2->v[ip2] += vtemp2;
              // s2->x[ip2] += update->dt * vtemp2;
            }
          }
        }
      }
    }
  }
}

// void FixContactHertz::post_advance_particles() {
//   // cout << "In FixContactHertz::initial_integrate()\n";

//   // Go through all the particles in the group and set b to the right value:
//   Vector3d f, dx;

//   Solid *s1, *s2;
//   Vector3d ftot, ftot_reduced, vtemp1, vtemp2;

//   double Rp, Rp1, Rp2, r, p, fmag, Estar, max_cellsize;

//   ftot = Vector3d();

//   s1 = domain->solids[solid1];
//   s2 = domain->solids[solid2];

//   Estar = 1.0 / ((1 - s1->mat->nu * s1->mat->nu) / s1->mat->E +
//                  (1 - s2->mat->nu * s2->mat->nu) / s2->mat->E);

//   max_cellsize = MAX(s1->grid->cellsize, s2->grid->cellsize);

//   if (domain->dimension == 2) {
//     for (int ip1 = 0; ip1 < s1->np_local; ip1++) {
//       for (int ip2 = 0; ip2 < s2->np_local; ip2++) {
//         dx = s2->x[ip2] - s1->x[ip1];

//     // Extremely gross screening:
//     if ((dx[0] < max_cellsize) && (dx[1] < max_cellsize) &&
//         (dx[2] < max_cellsize) && (dx[0] > -max_cellsize) &&
//             (dx[1] > -max_cellsize) && (dx[2] > -max_cellsize)) {
//       Rp1 = 0.5 * pow(s1->vol[ip1], 0.333333333);
//       Rp2 = 0.5 * pow(s2->vol[ip2], 0.333333333);
//       Rp = Rp1 + Rp2;


//       // Gross screening:
//       if ((dx[0] < Rp) && (dx[1] < Rp) && (dx[2] < Rp) && (dx[0] > -Rp) &&
//           (dx[1] > -Rp) && (dx[2] > -Rp)) {

//         r = dx.norm();

//         // Finer screening:
//         if (r < Rp) {

//               p = Rp - r; // penetration

//               fmag = 0.25 * M_PI * Estar *
//                      sqrt(Rp1 * Rp2 / (Rp1 + Rp2) * p * p * p);

//               f = fmag * dx / r;
//               ftot += f;
//           //s1->mbp[ip1] -= f;
//           //s2->mbp[ip2] += f;
//               vtemp1 = -update->dt * f / s1->mass[ip1];
//               s1->v[ip1] += vtemp1;
//               s1->x[ip1] += update->dt * vtemp1;

//               vtemp2 = update->dt * f / s2->mass[ip2];
//               s2->v[ip2] += vtemp2;
//               s2->x[ip2] += update->dt * vtemp2;
//             }
//           }
//         }
//       }
//     }
//   }
//   else if (domain->dimension == 3) {
//     for (int ip1 = 0; ip1 < s1->np_local; ip1++) {
//       for (int ip2 = 0; ip2 < s2->np_local; ip2++) {
//         dx = s2->x[ip2] - s1->x[ip1];

//     // Extremely gross screening:
//     if ((dx[0] < max_cellsize) && (dx[1] < max_cellsize) &&
//         (dx[2] < max_cellsize) && (dx[0] > -max_cellsize) &&
//             (dx[1] > -max_cellsize) && (dx[2] > -max_cellsize)) {
//       Rp1 = 0.5 * pow(s1->vol[ip1], 0.333333333);
//       Rp2 = 0.5 * pow(s2->vol[ip2], 0.333333333);
//       Rp = Rp1 + Rp2;


//       // Gross screening:
//       if ((dx[0] < Rp) && (dx[1] < Rp) && (dx[2] < Rp) && (dx[0] > -Rp) &&
//           (dx[1] > -Rp) && (dx[2] > -Rp)) {

//         r = dx.norm();

//         // Finer screening:
//         if (r < Rp) {
//               p = Rp - r; // penetration

//               fmag = four_thirds * Estar *
//                      sqrt(Rp1 * Rp2 / (Rp1 + Rp2) * p * p * p);

//               f = fmag * dx / r;
//           s1->mbp[ip1] -= f;
//           s2->mbp[ip2] += f;
//               ftot += f;
//               // vtemp1 = -update->dt * f / s1->mass[ip1];
//               // s1->v[ip1] += vtemp1;
//               // s1->x[ip1] += update->dt * vtemp1;

//               // vtemp2 = update->dt * f / s2->mass[ip2];
//               // s2->v[ip2] += vtemp2;
//               // s2->x[ip2] += update->dt * vtemp2;
//             }
//           }
//         }
//       }
//     }
//   }

//   // Reduce ftot:
//   MPI_Allreduce(ftot.elements, ftot_reduced.elements, 3, MPI_DOUBLE, MPI_SUM,
//                 universe->uworld);

//   (*input->vars)[id + "_x"] = Var(id + "_x", ftot_reduced[0]);
//   (*input->vars)[id + "_y"] = Var(id + "_y", ftot_reduced[1]);
//   (*input->vars)[id + "_z"] = Var(id + "_z", ftot_reduced[2]);
// }

void FixContactHertz::write_restart(ofstream *of) {
  of->write(reinterpret_cast<const char *>(&solid1), sizeof(int));
  of->write(reinterpret_cast<const char *>(&solid2), sizeof(int));
}

void FixContactHertz::read_restart(ifstream *ifr) {
  ifr->read(reinterpret_cast<char *>(&solid1), sizeof(int));
  ifr->read(reinterpret_cast<char *>(&solid2), sizeof(int));
}
