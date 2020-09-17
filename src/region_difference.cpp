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

#include "region_difference.h"
#include "domain.h"
#include "error.h"
#include "input.h"
#include "var.h"
#include <iostream>

using namespace std;

#define BIG 1.0e20

Difference::Difference(MPM *mpm, vector<string> args) : Region(mpm, args) {
  cout << "Initiate Difference" << endl;

  if (args.size() < Nargs) {
    error->all(FLERR,
               "Error: region_difference command not enough arguments\n" +
                   usage);
  }

  if (args.size() > Nargs) {
    error->all(FLERR,
               "Error: region_difference command too many arguments\n" +
                   usage);
  }

  vector<double> lim;
  int ir;

  // Look for the regionID:
  ir = domain->find_region(args[2]);

  if (ir == -1) {
    error->all(FLERR, "Error: region " + args[2] + " does not exist.\n");
  } else {
    lim.clear();
    iregions.push_back(ir);
    lim = domain->regions[ir]->limits();

    xlo = lim[0];
    xhi = lim[1];
    ylo = lim[2];
    yhi = lim[3];
    zlo = lim[4];
    zhi = lim[5];
  }

  ir = domain->find_region(args[3]);

  if (ir == -1) {
    error->all(FLERR, "Error: region " + args[3] + " does not exist.\n");
  } else {
    iregions.push_back(ir);
  }
}

Difference::~Difference()
{

}

/* ----------------------------------------------------------------------
   inside = 1 if x,y,z is inside or on surface
   inside = 0 if x,y,z is outside and not on surface
------------------------------------------------------------------------- */

int Difference::inside(double x, double y, double z) {
  if (domain->regions[iregions[0]]->inside(x, y, z) == 1 &&
      domain->regions[iregions[1]]->inside(x, y, z) == 0) {
    return 1;
  } else {
    return 0;
  }
}

/* ----------------------------------------------------------------------
   return a vector that contains the limits of the box
------------------------------------------------------------------------- */
vector<double> Difference::limits(){
  vector<double> lim;
  lim.push_back(xlo);
  lim.push_back(xhi);
  lim.push_back(ylo);
  lim.push_back(yhi);
  lim.push_back(zlo);
  lim.push_back(zhi);
  return lim;
}
