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

#include <fix_contact_min_penetration_plane.h>
#include <domain.h>
#include <error.h>
#include <group.h>
#include <input.h>
#include <solid.h>
#include <universe.h>
#include <update.h>

using namespace std;
using namespace FixConst;

FixContactMinPenetrationPlane::FixContactMinPenetrationPlane(MPM *mpm, vector<string> args)
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

    D = 0;
    xq = Vector3d();
    n = Vector3d();
    mu = 0;
    return;
  }

  if (args.size() < Nargs) {
    error->all(FLERR, "Error: not enough arguments.\n" + usage);
  }

  if (universe->me == 0) {
    cout << "Creating new fix FixContactMinPenetrationPlane with ID: " << args[0] << endl;
  }
  id = args[0];
  requires_ghost_particles = false;

  xq[0] = input->parsev(args[3]);
  xq[1] = input->parsev(args[4]);
  xq[2] = input->parsev(args[5]);

  n[0] = input->parsev(args[6]);
  n[1] = input->parsev(args[7]);
  n[2] = input->parsev(args[8]);

  mu = input->parsev(args[9]);

  // Normalize:
  n /= n.norm();
  
  D = -n[0] * xq[0] - n[1] * xq[1] - n[2] * xq[2];
}

void FixContactMinPenetrationPlane::prepare()
{
  ftot = Vector3d();
}

void FixContactMinPenetrationPlane::reduce()
{
  Vector3d ftot_reduced;

  // Reduce ftot:
  MPI_Allreduce(ftot.elements, ftot_reduced.elements, 3, MPI_FLOAT, MPI_SUM,
                universe->uworld);

  (*input->vars)[id + "_x"] = Var(id + "_x", ftot_reduced[0]);
  (*input->vars)[id + "_y"] = Var(id + "_y", ftot_reduced[1]);
  (*input->vars)[id + "_z"] = Var(id + "_z", ftot_reduced[2]);
}

void FixContactMinPenetrationPlane::initial_integrate(Solid &solid)
{
  // cout << "In FixContactMinPenetrationPlane::initial_integrate()\n";

  // Go through all the particles in the group and set b to the right value:

  for (int ip = 0; ip < solid.np_local; ip++)
  {
    // Extremely gross screening:
    float d = n.dot(solid.x[ip] - xq);

    if (d >= solid.grid->cellsize)
      continue;

    float Rp;
    if (domain->dimension == 2)
    {
      Rp = sqrt(solid.vol[ip]);
    }
    else
    {
      Rp = cbrt(solid.vol[ip]);
    }
  
    // Fine screening:
    float p = Rp/2 - d;

    if (p < 0)
      continue;

    float fnorm = solid.mass[ip]*p/update->dt/update->dt;
    Vector3d f = fnorm*n;

    if (mu)
    {
      Vector3d vt = solid.v[ip] - n.dot(solid.v[ip])*n;
      float vtnorm = vt.norm();
      if (vtnorm)
      {
        vt /= vtnorm;
        f -= mu*fnorm*vt;
      }
    }

    solid.mbp[ip] += f;
    ftot += f;
  }
}

void FixContactMinPenetrationPlane::write_restart(ofstream *of)
{
  of->write(reinterpret_cast<const char *>(&D), sizeof(float));
  of->write(reinterpret_cast<const char *>(&mu), sizeof(float));
  of->write(reinterpret_cast<const char *>(&xq), sizeof(Vector3d));
  of->write(reinterpret_cast<const char *>(&n), sizeof(Vector3d));
}

void FixContactMinPenetrationPlane::read_restart(ifstream *ifr)
{
  ifr->read(reinterpret_cast<char *>(&D), sizeof(float));
  ifr->read(reinterpret_cast<char *>(&mu), sizeof(float));
  ifr->read(reinterpret_cast<char *>(&xq), sizeof(Vector3d));
  ifr->read(reinterpret_cast<char *>(&n), sizeof(Vector3d));
}
