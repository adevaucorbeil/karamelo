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

#include <fix_contact_pinball.h>
#include <domain.h>
#include <error.h>
#include <group.h>
#include <input.h>
#include <solid.h>
#include <universe.h>
#include <update.h>

using namespace std;
using namespace FixConst;

FixContactPinball::FixContactPinball(MPM *mpm, vector<string> args)
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

    K = 0;
    return;
  }

  if (args.size() < Nargs) {
    error->all(FLERR, "Error: not enough arguments.\n" + usage);
  }

  if (universe->me == 0) {
    cout << "Creating new fix FixContactPinball with ID: " << args[0] << endl;
  }
  id = args[0];
  requires_ghost_particles = true;

  K = input->parsev(args[4]);
}

void FixContactPinball::prepare()
{
  ftot = Vector3d();
}

void FixContactPinball::reduce()
{
  Vector3d ftot_reduced;

  // Reduce ftot:
  MPI_Allreduce(ftot.elements, ftot_reduced.elements, 3, MPI_DOUBLE, MPI_SUM,
                universe->uworld);

  (*input->vars)[id + "_x"] = Var(id + "_x", ftot_reduced[0]);
  (*input->vars)[id + "_y"] = Var(id + "_y", ftot_reduced[1]);
  (*input->vars)[id + "_z"] = Var(id + "_z", ftot_reduced[2]);
}

void FixContactPinball::initial_integrate(Solid &solid)
{
  // cout << "In FixContactPinball::initial_integrate()\n";

  // Go through all the particles in the group and set b to the right value:
  for (int ip = 0; ip < solid.np_local; ip++)
    for (Solid *solid1: domain->solids)
    {
      if (solid1 <= &solid)
        continue;
    
      double Gstar = 1/((1 - solid.  mat->nu*solid.  mat->nu)/solid.  mat->G +
                        (1 - solid1->mat->nu*solid1->mat->nu)/solid1->mat->G);

      double max_cellsize = MAX(solid.grid->cellsize, solid1->grid->cellsize);

      for (int ip1 = 0; ip1 < solid1->np_local; ip1++)
      {
        const Vector3d &dx = solid.x[ip1] - solid.x[ip];
      
        // Extremely gross screening:
        bool outside = false;

        for (int i = 0; i < domain->dimension; i++)
          if (outside = abs(dx[i]) > max_cellsize)
            break;

        if (outside)
          continue;

        double Rp0, Rp1;
        if (domain->dimension == 2)
        {
          if (domain->axisymmetric)
          {
            Rp0 = sqrt(solid.  vol[ip ]/solid.x[ip ][0])/2;
            Rp1 = sqrt(solid1->vol[ip1]/solid.x[ip1][0])/2;
          }
          else
          {
            Rp0 = sqrt(solid.  vol[ip ])/2;
            Rp1 = sqrt(solid1->vol[ip1])/2;
          }
        }
        else
        {
          Rp0 = cbrt(solid.  vol[ip ])/2;
          Rp1 = cbrt(solid1->vol[ip1])/2;
        }
        double Rp = Rp0 + Rp1;
      
        // Gross screening:
        for (int i = 0; i < domain->dimension; i++)
          if (outside = abs(dx[i]) > Rp)
            break;

        if (outside)
          continue;
      
        double r = dx.norm();

        // Finer screening:
        if (r >= Rp)
          continue;

        double p = Rp - r; // penetration

        const Vector3d &dv = solid1->v[ip1] - solid.v[ip];
        double pdot = -dx.dot(dv)/r; // penetration speed

        if (pdot < 0)
          continue;

        double Rp03 = Rp0*Rp0*Rp0;
        double Rp13 = Rp1*Rp1*Rp1;

        const Vector3d &f = MIN(solid.rho[ip]*Rp03 * solid1->rho[ip1]*Rp13*pdot/
                               (solid.rho[ip]*Rp03 + solid1->rho[ip1]*Rp13)/update->dt,
                                Gstar*sqrt(Rp0*Rp1/(Rp0 + Rp1)*p)*p)*dx/r;
        ftot += f;
        solid.  mbp[ip ] -= f;
        solid1->mbp[ip1] += f;
      }
    }
}

void FixContactPinball::write_restart(ofstream *of)
{
  of->write(reinterpret_cast<const char *>(&K), sizeof(double));
}

void FixContactPinball::read_restart(ifstream *ifr)
{
  ifr->read(reinterpret_cast<char *>(&K), sizeof(double));
}
