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

#include <region_sphere.h>
#include <domain.h>
#include <error.h>
#include <input.h>
#include <math_special.h>
#include <universe.h>
#include <update.h>
#include <var.h>
#include <iostream>

using namespace std;
using namespace MathSpecial;

#define BIG 1.0e20

Sphere::Sphere(MPM *mpm, vector<string> args) : Region(mpm, args) {
  if (universe->me == 0)
    cout << "Initiate Sphere" << endl;

  c1 = c2 = c3 = 0;
  R = 0;

  if (args.size() < 3) {
    error->all(FLERR, "Error: not enough arguments.\n");
  }

  if (args[2].compare("restart") ==
      0) { // If the keyword restart, we are expecting to have read_restart()
           // launched right after.
    RSq = xlo = xhi = ylo = yhi = zlo = zhi = 0;
    return;
  }

  if (args.size() < Nargs[domain->dimension - 1]) {
    error->all(FLERR, "Error: region command not enough arguments.\n" +
                          usage[domain->dimension - 1]);
  }

  if (args.size() > Nargs[domain->dimension - 1]) {
    error->all(FLERR, "Error: region command too many arguments.\n" +
                          usage[domain->dimension - 1]);
  }

  int iargs = 2;
  c1 = input->parsev(args[iargs++]);
  if (domain->dimension >= 2)
    c2 = input->parsev(args[iargs++]);
  if (domain->dimension == 3)
    c3 = input->parsev(args[iargs++]);
  R = input->parsev(args[iargs++]);

  // Chech if R is negative:
  if (R < 0) {
    error->all(FLERR,
               "Error: R cannot be negative, set to: " + to_string(R) + "\n.");
  }

  RSq = R * R;
  if (universe->me == 0)
    cout << "c1, c2, c3, R = " << c1 << "\t" << c2 << "\t" << c3 << "\t" << R
         << endl;

  xlo = c1 - R;
  xhi = c1 + R;
  ylo = c2 - R;
  yhi = c2 + R;
  zlo = c3 - R;
  zhi = c3 + R;

  if (update->method_type == "tlmpm") {
    if (domain->boxlo[0] > xlo)
      domain->boxlo[0] = xlo;
    if (domain->boxlo[1] > ylo)
      domain->boxlo[1] = ylo;
    if (domain->dimension == 3)
      if (domain->boxlo[2] > zlo)
        domain->boxlo[2] = zlo;

    if (domain->boxhi[0] < xhi)
      domain->boxhi[0] = xhi;
    if (domain->boxhi[1] < yhi)
      domain->boxhi[1] = yhi;
    if (domain->dimension == 3)
      if (domain->boxhi[2] < zhi)
        domain->boxhi[2] = zhi;
  } else {
    if (domain->boxlo[0] > xlo)
      xlo = domain->boxlo[0];
    if (domain->boxlo[1] > ylo)
      ylo = domain->boxlo[1];
    if (domain->dimension == 3)
      if (domain->boxlo[2] > zlo)
        zlo = domain->boxlo[2];
  }
}

Sphere::~Sphere()
{

}

/* ----------------------------------------------------------------------
   inside = 1 if x,y,z is inside or on surface
   inside = 0 if x,y,z is outside and not on surface
------------------------------------------------------------------------- */

int Sphere::inside(float x, float y, float z)
{
  //cout << "Check if point (" << x << ", " << y << ", " << z << ") is inside the region" << endl;
  float dSq;
  dSq = square(x - c1) + square(y - c2) + square(z - c3);
  if (dSq <= RSq) {
    return 1;
  }
  return 0;
}

/* ----------------------------------------------------------------------
   return a vector that contains the limits of the box
------------------------------------------------------------------------- */
vector<float> Sphere::limits(){
  vector<float> lim;
  lim.push_back(xlo);
  lim.push_back(xhi);
  lim.push_back(ylo);
  lim.push_back(yhi);
  lim.push_back(zlo);
  lim.push_back(zhi);
  return lim;
}

void Sphere::write_restart(ofstream *of) {
  of->write(reinterpret_cast<const char *>(&c1), sizeof(float));
  of->write(reinterpret_cast<const char *>(&c2), sizeof(float));
  of->write(reinterpret_cast<const char *>(&c3), sizeof(float));
  of->write(reinterpret_cast<const char *>(&R), sizeof(float));
  of->write(reinterpret_cast<const char *>(&xlo), sizeof(float));
  of->write(reinterpret_cast<const char *>(&xhi), sizeof(float));
  of->write(reinterpret_cast<const char *>(&ylo), sizeof(float));
  of->write(reinterpret_cast<const char *>(&yhi), sizeof(float));
  of->write(reinterpret_cast<const char *>(&zlo), sizeof(float));
  of->write(reinterpret_cast<const char *>(&zhi), sizeof(float));
}



void Sphere::read_restart(ifstream *ifr) {
  // cout << "Restart Sphere" << endl;
  ifr->read(reinterpret_cast<char *>(&c1), sizeof(float));
  ifr->read(reinterpret_cast<char *>(&c2), sizeof(float));
  ifr->read(reinterpret_cast<char *>(&c3), sizeof(float));
  ifr->read(reinterpret_cast<char *>(&R), sizeof(float));
  RSq = R * R;
  ifr->read(reinterpret_cast<char *>(&xlo), sizeof(float));
  ifr->read(reinterpret_cast<char *>(&xhi), sizeof(float));
  ifr->read(reinterpret_cast<char *>(&ylo), sizeof(float));
  ifr->read(reinterpret_cast<char *>(&yhi), sizeof(float));
  ifr->read(reinterpret_cast<char *>(&zlo), sizeof(float));
  ifr->read(reinterpret_cast<char *>(&zhi), sizeof(float));
}
