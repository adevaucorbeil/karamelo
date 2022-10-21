/* ----------------------------------------------------------------------
 *
 *                    ***       Karamelo       ***
 *               Parallel Material Point Method Simulator
 *
 * Coprright (2019) Alban de Vaucorbeil, alban.devaucorbeil@monash.edu
 * Materials Science and Engineering, Monash University
 * Clayton VIC 3800, Australia

 * This software is distributed under the GNU General Public License.
 *
 * ----------------------------------------------------------------------- */

#include "fix_contact.h"
#include "domain.h"
#include "error.h"
#include "group.h"
#include "input.h"
#include "method.h"
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

// is absolute of x less than absolute of r?
#define ABS_LESS(x, r) (x < r) && (x > -r)
#define ABS_LESS_XY(x, r) ABS_LESS(x[0], r) && ABS_LESS(x[0], r)
#define ABS_LESS_XYZ(x, r) ABS_LESS_XY(x, r) && ABS_LESS(x[2], r)

FixContact::FixContact(MPM *mpm, vector<string> args)
    : Fix(mpm, args)
{
  if (args.size() < 3)
  {
    error->all(FLERR, "Error: not enough arguments.\n");
  }

  if (args[2].compare("restart") == 0)
  { // If the keyword restart, we are expecting to have read_restart()
    // launched right after.
    igroup = stoi(args[3]);
    if (igroup == -1 && universe->me == 0)
    {
      cout << "Could not find group number " << args[3] << endl;
    }
    groupbit = group->bitmask[igroup];

    solid1 = solid2 = -1;
    return;
  }

  if (args.size() < this->NARGS)
  {
    error->all(FLERR, "Error: not enough arguments.\n" + this->USAGE);
  }

  solid1 = domain->find_solid(args[2]);
  if (solid1 < 0)
  {
    error->all(FLERR, "Error: solid " + args[2] + " unknown.\n");
  }

  solid2 = domain->find_solid(args[3]);
  if (solid1 < 0)
  {
    error->all(FLERR, "Error: solid " + args[3] + " unknown.\n");
  }

  if (universe->me == 0)
  {
    cout << "Creating new fix FixContact " << args[1] << " with ID: " << args[0] << endl;
  }
  id = args[0];
  requires_ghost_particles = true;
}

void FixContact::init() {
  Solid * s1, * s2;
  s1 = domain->solids[solid1];
  s2 = domain->solids[solid2];
  double max_cellsize;
  max_cellsize = MAX(s1->grid->cellsize, s2->grid->cellsize);
  // cout << "init solids " << solid1 << " and " << solid2 << endl;
  s1->surfmask_init(max_cellsize);
  s2->surfmask_init(max_cellsize);
  // cout << "init surfmask of solids " << solid1 << " and " << solid2 << endl;
}

void FixContact::setup() {}

void FixContact::setmask()
{
  mask = 0;
  mask |= INITIAL_INTEGRATE;
  // mask |= POST_ADVANCE_PARTICLES;
}

void FixContact::initial_integrate()
{
  // cout << "In FixContact::initial_integrate()\n";

  // Go through all the particles in the group and set b to the right value:
  Solid *s1, *s2;
  Eigen::Vector3d dx, ftot, ftot_reduced, vtemp1, vtemp2, dv, vt;

  double Rp, Rp1, Rp2, r, max_cellsize;

  ftot.setZero();

  s1 = domain->solids[solid1];
  s2 = domain->solids[solid2];

  int root = 0;
  // it could be better to gather every particle to that process in which domain it is at the moment
  // instead of gathering all on the same
  // s1->gather_particles(root);
  // s2->gather_particles(root);
  // in the future, when there could be more than 2 solids in contact, it will be necessary to update the max cellsize when 
  // a new fix contact is added
  max_cellsize = MAX(s1->get_surfcellsize(), s2->get_surfcellsize());
  // cout << "solid: " << solid1 << endl;
  s1->compute_surface_particles();
  // cout << "solid: " << solid2 << endl;
  s2->compute_surface_particles();

  if (domain->dimension == 2)
  {
    for (int ip1 = 0; ip1 < s1->np_local; ip1++)
    {
      if (!s1->is_surf[ip1])
        continue;

      for (int ip2 = 0; ip2 < s2->np_local; ip2++)
      {
        if (!s2->is_surf[ip2])
          continue;

        dx = s2->x[ip2] - s1->x[ip1];

        // Extremely gross screening:
        if (ABS_LESS_XY(dx, max_cellsize))
        {
          if (domain->axisymmetric)
          {
            Rp1 = 0.5 * sqrt(s1->vol[ip1] / s1->x[ip1][0]);
            Rp2 = 0.5 * sqrt(s2->vol[ip2] / s2->x[ip2][0]);
          }
          else
          {
            Rp1 = 0.5 * sqrt(s1->vol[ip1]);
            Rp2 = 0.5 * sqrt(s2->vol[ip2]);
          }
          Rp = Rp1 + Rp2;

          // Gross screening:
          if (ABS_LESS_XY(dx, Rp))
          {

            r = dx.norm();

            // Finer screening:
            if (r < Rp)
            {
              force_increment(dx, ftot, s1, s2, ip1, ip2, r, Rp1, Rp2);
            }
          }
        }
      }
    }
  }
  else if (domain->dimension == 3)
  {
    for (int ip1 = 0; ip1 < s1->np_local; ip1++)
    {
      for (int ip2 = 0; ip2 < s2->np_local; ip2++)
      {
        dx = s2->x[ip2] - s1->x[ip1];

        // Extremely gross screening:
        if (ABS_LESS_XYZ(dx, max_cellsize))
        {
          Rp1 = 0.5 * cbrt(s1->vol[ip1]);
          Rp2 = 0.5 * cbrt(s2->vol[ip2]);
          Rp = Rp1 + Rp2;

          // Gross screening:
          if (ABS_LESS_XYZ(dx, Rp))
          {

            r = dx.norm();
            // Finer screening:
            if (r < Rp)
            {
              force_increment(dx, ftot, s1, s2, ip1, ip2, r, Rp1, Rp2);
            }
          }
        }
      }
    }
  }

  // s1->scatter_particles(root);
  // s2->scatter_particles(root);

  // Reduce ftot:
  MPI_Allreduce(ftot.data(), ftot_reduced.data(), 3, MPI_DOUBLE, MPI_SUM,
                universe->uworld);

  (*input->vars)[id + "_x"] = Var(id + "_x", ftot_reduced[0]);
  (*input->vars)[id + "_y"] = Var(id + "_y", ftot_reduced[1]);
  (*input->vars)[id + "_z"] = Var(id + "_z", ftot_reduced[2]);
}

// void FixContact::post_advance_particles() {
//   // cout << "In FixContact::initial_integrate()\n";

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

void FixContact::write_restart(ofstream *of)
{
  of->write(reinterpret_cast<const char *>(&solid1), sizeof(int));
  of->write(reinterpret_cast<const char *>(&solid2), sizeof(int));
}

void FixContact::read_restart(ifstream *ifr)
{
  ifr->read(reinterpret_cast<char *>(&solid1), sizeof(int));
  ifr->read(reinterpret_cast<char *>(&solid2), sizeof(int));
}
