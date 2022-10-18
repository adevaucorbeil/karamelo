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

#include "fix_contact_min_penetration.h"
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

FixContactMinPenetration::FixContactMinPenetration(MPM *mpm, vector<string> args)
    : FixContact(mpm, args)
{
  Solid *s1, *s2;
  s1 = domain->solids[solid1];
  s2 = domain->solids[solid2];

  alpha = s1->mat->kappa / (s1->mat->kappa + s2->mat->kappa);
  mu = input->parsev(args[4]);
}

FixContactMinPenetration::~FixContactMinPenetration(){};

void FixContactMinPenetration::force_increment(
    Eigen::Vector3d &dx, Eigen::Vector3d &ftot,
    Solid *s1, Solid *s2,
    const int ip1, const int ip2,
    const double r, const double Rp1, const double Rp2)
{
  double fmag, ffric, inv_r, vtnorm, gamma;
  Eigen::Vector3d dv, vt, f;
  inv_r = 1.0 / r;
  fmag = s1->mass[ip1] * s2->mass[ip2] / ((s1->mass[ip1] + s2->mass[ip2]) * update->dt * update->dt) * (1 - (Rp1 + Rp2) * inv_r);
  f = fmag * dx;

  if (mu != 0)
  {
    dv = s2->v[ip2] - s1->v[ip1];
    vt = dv - dv.dot(dx) * inv_r * inv_r * dx;
    vtnorm = vt.norm();
    if (vtnorm != 0)
    {
      vt /= vtnorm;
      ffric = mu * fmag * r;
      f -= ffric * vt;
      if (update->method->temp)
      {
        gamma = ffric * vtnorm * update->dt;
        s1->gamma[ip1] += alpha * s1->vol0[ip1] * s1->mat->invcp * gamma;
        s2->gamma[ip2] += (1.0 - alpha) * s2->vol0[ip2] * s2->mat->invcp * gamma;
      }
    }
  }

  s1->mbp[ip1] += f;
  s2->mbp[ip2] -= f;
  ftot += f;

  // if ((s1->ptag[ip1] == 12862 && s2->ptag[ip2] == 12902) ||
  //     (s1->ptag[ip1] == 12835 && s2->ptag[ip2] == 12875)) {
  //   cout << "dx[" << s1->ptag[ip1] << "," << s2->ptag[ip2] << "]= ["
  //        << dx[0] << "\t" << dx[1] << "\t" << dx[2] << "]\n";
  //   cout << "f[" << s1->ptag[ip1] << "," << s2->ptag[ip2] << "]= ["
  //        << f[0] << "\t" << f[1] << "\t" << f[2] << "]\n";
  // }
}

void FixContactMinPenetration::initial_integrate()
{
  // cout << "In FixContactMinPenetration::initial_integrate()\n";

  // Go through all the particles in the group and set b to the right value:
  Eigen::Vector3d f, dx;

  Solid *s1, *s2;

  s1 = domain->solids[solid1];
  s2 = domain->solids[solid2];

  FixContact::initial_integrate();
}

void FixContactMinPenetration::write_restart(ofstream *of)
{
  FixContact::write_restart(of);
  of->write(reinterpret_cast<const char *>(&mu), sizeof(double));
}

void FixContactMinPenetration::read_restart(ifstream *ifr)
{
  FixContact::read_restart(ifr);
  ifr->read(reinterpret_cast<char *>(&mu), sizeof(double));
}
