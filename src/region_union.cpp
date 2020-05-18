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

#include "region_union.h"
#include "domain.h"
#include "error.h"
#include "input.h"
#include "var.h"
#include <iostream>

using namespace std;

#define BIG 1.0e20

Union::Union(MPM *mpm, vector<string> args) : Region(mpm, args)
{
  cout << "Initiate Union" << endl;

  if (args.size() < Nargs) {
    error->all(FLERR,
               "Error: region_union command not enough arguments\n" + usage);
  }

  vector<double> lim;

  xlo = BIG;
  xhi = -BIG;
  ylo = BIG;
  yhi = -BIG;
  zlo = BIG;
  zhi = -BIG;

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

      if (xlo > lim[0]) xlo = lim[0];
      if (xhi < lim[1]) xhi = lim[1];
      if (ylo > lim[2]) ylo = lim[2];
      if (yhi < lim[3]) yhi = lim[3];
      if (zlo > lim[4]) zlo = lim[4];
      if (zhi < lim[5]) zhi = lim[5];
    }
  }
}


Union::~Union()
{

}

/* ----------------------------------------------------------------------
   inside = 1 if x,y,z is inside or on surface
   inside = 0 if x,y,z is outside and not on surface
------------------------------------------------------------------------- */

int Union::inside(double x, double y, double z)
{
  for(int i; i < iregions.size(); i++) {
    if (domain->regions[iregions[i]]->inside(x, y, z) == 1)
      return 1;
  }
  return 0;
}

/* ----------------------------------------------------------------------
   return a vector that contains the limits of the box
------------------------------------------------------------------------- */
vector<double> Union::limits(){
  vector<double> lim;
  lim.push_back(xlo);
  lim.push_back(xhi);
  lim.push_back(ylo);
  lim.push_back(yhi);
  lim.push_back(zlo);
  lim.push_back(zhi);
  return lim;
}
