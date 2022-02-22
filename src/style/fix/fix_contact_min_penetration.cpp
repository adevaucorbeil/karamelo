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

#include <fix_contact_min_penetration.h>
#include <domain.h>
#include <error.h>
#include <group.h>
#include <input.h>
#include <method.h>
#include <solid.h>
#include <universe.h>
#include <update.h>

using namespace std;
using namespace FixConst;

FixContactMinPenetration::FixContactMinPenetration(MPM *mpm, vector<string> args):
  Fix(mpm, args, INITIAL_INTEGRATE)
{
  if (args.size() < 3)
    error->all(FLERR, "Error: not enough arguments.\n");

  if (args[2] == "restart")
  {
    // If the keyword restart, we are expecting to have read_restart() launched right after.
    igroup = stoi(args[3]);
    if (igroup == -1 && !universe->me)
      cout << "Could not find group number " << args[3] << endl;

    groupbit = group->bitmask[igroup];

    mu = 0;
    return;
  }

  if (args.size() < Nargs)
    error->all(FLERR, "Error: not enough arguments.\n" + usage);

  if (!universe->me)
    cout << "Creating new fix FixContactMinPenetration with ID: " << args[0] << endl;

  id = args[0];
  requires_ghost_particles = true;

  mu = input->parsev(args[4]);
}

void FixContactMinPenetration::prepare()
{
  ftot = Vector3d();
}

void FixContactMinPenetration::reduce()
{
  Vector3d ftot_reduced;

  // Reduce ftot:
  MPI_Allreduce(ftot.elements, ftot_reduced.elements, 3, MPI_DOUBLE, MPI_SUM,
                universe->uworld);

  (*input->vars)[id + "_x"] = Var(id + "_x", ftot_reduced[0]);
  (*input->vars)[id + "_y"] = Var(id + "_y", ftot_reduced[1]);
  (*input->vars)[id + "_z"] = Var(id + "_z", ftot_reduced[2]);
}

void FixContactMinPenetration::initial_integrate(Solid &solid, int ip)
{
  // cout << "In FixContactMinPenetration::initial_integrate()\n";

  // Go through all the particles in the group and set b to the right value:
  for (Solid *solid1: domain->solids)
  {
    if (solid1 <= &solid)
      continue;

    double Estar = 1/((1 - solid.  mat->nu*solid.  mat->nu)/solid.  mat->E +
                      (1 - solid1->mat->nu*solid1->mat->nu)/solid1->mat->E);

    double alpha = solid.mat->kappa/(solid.mat->kappa + solid1->mat->kappa);

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

      double Rp;
      if (domain->dimension == 2)
      {
        if (domain->axisymmetric)
          Rp = (sqrt(solid.  vol.at(ip )/solid.x[ip ][0]) +
                sqrt(solid1->vol.at(ip1)/solid.x[ip1][0]))/2;
        else
          Rp = (sqrt(solid.vol.at(ip)) + sqrt(solid1->vol.at(ip1)))/2;
      }
      else
        Rp = (cbrt(solid.vol.at(ip)) + cbrt(solid1->vol.at(ip1)))/2;
      
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

      double fmag = solid.mass.at(ip) * solid1->mass.at(ip1)/
                  ((solid.mass.at(ip) + solid1->mass.at(ip1))*
                   update->dt*update->dt)*(1 - Rp/r);
      Vector3d f = fmag*dx;

      if (mu)
      {
        const Vector3d &dv = solid1->v.at(ip1) - solid.v.at(ip);
        Vector3d vt = dv - dv.dot(dx)/r/r*dx;
        double vtnorm = vt.norm();
        if (vtnorm != 0)
        {
          vt /= vtnorm;
          double ffric = mu*fmag*r;
          f -= ffric*vt;

          if (update->method->temp)
          {
            double gamma = ffric*vtnorm*update->dt;
            solid.gamma.at(ip) += alpha*solid.vol0.at(ip)*solid.mat->invcp*gamma;
            solid1->gamma.at(ip1) += (1 - alpha)*solid1->vol0.at(ip1)*solid1->mat->invcp*gamma;
          }
        }
      }

      solid.  mbp.at(ip ) += f;
      solid1->mbp.at(ip1) -= f;
      ftot += f;
    }
  }
}

void FixContactMinPenetration::write_restart(ofstream *of)
{
  of->write(reinterpret_cast<const char *>(&mu), sizeof(double));
}

void FixContactMinPenetration::read_restart(ifstream *ifr)
{
  ifr->read(reinterpret_cast<char *>(&mu), sizeof(double));
}
