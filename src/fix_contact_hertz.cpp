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

#include "fix_contact_hertz.h"
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

FixContactHertz::FixContactHertz(MPM *mpm, vector<string> args)
    : FixContact(mpm, args)
{
  Solid *s1, *s2;
  s1 = domain->solids[solid1];
  s2 = domain->solids[solid2];
  this->Estar = 1.0 / ((1 - s1->mat->nu * s1->mat->nu) / s1->mat->E +
                       (1 - s2->mat->nu * s2->mat->nu) / s2->mat->E);
}

FixContactHertz::~FixContactHertz(){};

void FixContactHertz::force_increment(
    Eigen::Vector3d &dx, Eigen::Vector3d &f, Eigen::Vector3d &ftot,
    Solid *s1, Solid *s2,
    const int ip1, const int ip2,
    const double r, const double Rp1, const double Rp2)
{

  double fmag, p;
  p = Rp1 + Rp2 - r; // penetration

  switch (domain->dimension)
  {
  case 2:
    fmag = 0.25 * M_PI * Estar * sqrt(Rp1 * Rp2 / (Rp1 + Rp2) * p * p * p);
    break;

  case 3:
    fmag = four_thirds * Estar * sqrt(Rp1 * Rp2 / (Rp1 + Rp2) * p * p * p);
    break;

  default:
    break;
  }

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

void FixContactHertz::initial_integrate()
{
  // cout << "In FixContactHertz::initial_integrate()\n";

  // Go through all the particles in the group and set b to the right value:
  Eigen::Vector3d f, dx;

  Solid *s1, *s2;

  s1 = domain->solids[solid1];
  s2 = domain->solids[solid2];

  double Estar = 1.0 / ((1 - s1->mat->nu * s1->mat->nu) / s1->mat->E +
                        (1 - s2->mat->nu * s2->mat->nu) / s2->mat->E);
  FixContact::initial_integrate();
}
