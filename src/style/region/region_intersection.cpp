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

#include <region_intersection.h>
#include <domain.h>
#include <error.h>
#include <input.h>
#include <universe.h>
#include <var.h>
#include <iostream>

using namespace std;

#define BIG 1.0e20

Intersection::Intersection(MPM *mpm, vector<string> args) : Region(mpm, args)
{
  if (universe->me == 0)
    cout << "Initiate Intersection" << endl;

  if (args.size() < 3) {
    error->all(FLERR, "Error: not enough arguments.\n");
  }

  if (args[2].compare("restart") ==
      0) { // If the keyword restart, we are expecting to have read_restart()
           // launched right after.
    xlo = xhi = ylo = yhi = zlo = zhi = 0;
    return;
  }

  if (args.size() < Nargs) {
    error->all(FLERR,
               "Error: region_intersection command not enough arguments\n" + usage);
  }

  vector<float> lim;

  xlo = -BIG;
  xhi = BIG;
  ylo = -BIG;
  yhi = BIG;
  zlo = -BIG;
  zhi = BIG;

  int ir;
  for (int i = 2; i < args.size(); i++) {
    // Look for the regionID:
    ir = domain->find_region(args[i]);

    if (ir == -1) {
      error->all(FLERR, "Error: region " + args[i] + " does not exist.\n");
    } else {
      lim.clear();
      iregions.push_back(ir);
      lim = domain->regions[ir]->limits();

      if (xlo < lim[0]) xlo = lim[0];
      if (xhi > lim[1]) xhi = lim[1];
      if (ylo < lim[2]) ylo = lim[2];
      if (yhi > lim[3]) yhi = lim[3];
      if (zlo < lim[4]) zlo = lim[4];
      if (zhi > lim[5]) zhi = lim[5];
    }
  }
}


Intersection::~Intersection()
{

}

/* ----------------------------------------------------------------------
   inside = 1 if x,y,z is inside or on surface
   inside = 0 if x,y,z is outside and not on surface
------------------------------------------------------------------------- */

int Intersection::inside(float x, float y, float z)
{
  for(int i = 0; i < iregions.size(); i++) {
    if (domain->regions[iregions[i]]->inside(x, y, z) == 0)
      return 0;
  }
  return 1;
}

/* ----------------------------------------------------------------------
   return a vector that contains the limits of the box
------------------------------------------------------------------------- */
vector<float> Intersection::limits(){
  vector<float> lim;
  lim.push_back(xlo);
  lim.push_back(xhi);
  lim.push_back(ylo);
  lim.push_back(yhi);
  lim.push_back(zlo);
  lim.push_back(zhi);
  return lim;
}

void Intersection::write_restart(ofstream *of) {
  size_t N = iregions.size();
  of->write(reinterpret_cast<const char *>(&N), sizeof(size_t));
  for (int i = 0; i < N; i++) {
    of->write(reinterpret_cast<const char *>(&iregions[i]), sizeof(int));
  }

  of->write(reinterpret_cast<const char *>(&xlo), sizeof(float));
  of->write(reinterpret_cast<const char *>(&xhi), sizeof(float));
  of->write(reinterpret_cast<const char *>(&ylo), sizeof(float));
  of->write(reinterpret_cast<const char *>(&yhi), sizeof(float));
  of->write(reinterpret_cast<const char *>(&zlo), sizeof(float));
  of->write(reinterpret_cast<const char *>(&zhi), sizeof(float));
}

void Intersection::read_restart(ifstream *ifr) {
  // cout << "Restart Intersection" << endl;

  size_t N;
  ifr->read(reinterpret_cast<char *>(&N), sizeof(size_t));
  iregions.resize(N);

  for (int i = 0; i < N; i++) {
  ifr->read(reinterpret_cast<char *>(&iregions[i]), sizeof(int));
  }

  ifr->read(reinterpret_cast<char *>(&xlo), sizeof(float));
  ifr->read(reinterpret_cast<char *>(&xhi), sizeof(float));
  ifr->read(reinterpret_cast<char *>(&ylo), sizeof(float));
  ifr->read(reinterpret_cast<char *>(&yhi), sizeof(float));
  ifr->read(reinterpret_cast<char *>(&zlo), sizeof(float));
  ifr->read(reinterpret_cast<char *>(&zhi), sizeof(float));
}
