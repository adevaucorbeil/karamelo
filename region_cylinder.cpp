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

#include "region_cylinder.h"
#include "domain.h"
#include "error.h"
#include "input.h"
#include "math_special.h"
#include "method.h"
#include "update.h"
#include "var.h"
#include <iostream>

using namespace std;
using namespace MathSpecial;

#define BIG 1.0e20

RegCylinder::RegCylinder(MPM *mpm, vector<string> args) : Region(mpm, args)
{
  cout << "Initiate RegCylinder" << endl;

  if (domain->dimension == 3) {
    if (args.size()<8) {
      error->all(FLERR, "Error: not enough arguments.\n");
    }

    options(&args, args.begin()+8);
  } else {
    if (args.size()<5) {
      error->all(FLERR, "Error: not enough arguments.\n");
    }    

    options(&args, args.begin()+5);
  }

  if (domain->dimension == 3) {
    if (args[2].compare("x") == 0) {
      axis = 'x';
    } else if (args[2].compare("y") == 0) {
      axis = 'y';
    } else if (args[2].compare("z") == 0) {
      axis = 'z';
    } else {
      error->all(FLERR, "Error: region cylinder axis not understood, expect x, y, or z, received " + args[2]+".\n");
    }

    c1 = input->parsev(args[3]);
    c2 = input->parsev(args[4]);
    R = input->parsev(args[5]);

    if (args[6].compare("INF") == 0 || args[6].compare("EDGE") == 0) {
      if (domain->regions.size() == 0) {
	error->all(FLERR, "Cannot use region INF or EDGE when box does not exist.\n");
      }
      lo = -BIG;
    } else lo = input->parsev(args[6]);

    if (args[7].compare("INF") == 0 || args[7].compare("EDGE") == 0) {
      if (domain->regions.size() == 0) {
	error->all(FLERR, "Cannot use region INF or EDGE when box does not exist.\n");
      }
      hi = BIG;
    } else hi = input->parsev(args[7]);

    cout << "lo hi = " << lo << "\t" << hi << endl;
  } else {
    axis = 'z';

    c1 = input->parsev(args[2]);
    c2 = input->parsev(args[3]);
    R = input->parsev(args[4]);
    lo = hi = 0;
  }

  RSq = R*R;
  cout << "axis, c1, c2, R = " << axis << "\t" << c1 << "\t" << c2 << "\t" << R << endl;

  // error check

  if (lo > hi) {
    error->all(FLERR, "Illegal region cylinder command: low is higher than high.\n");
  }

  if (axis=='x') {
    xlo = lo;
    xhi = hi;
    ylo = c1 - R;
    yhi = c1 + R;
    zlo = c2 - R;
    zhi = c2 + R;
  } else if (axis=='y') {
    xlo = c1 - R;
    xhi = c1 + R;
    ylo = lo;
    yhi = hi;
    zlo = c2 - R;
    zhi = c2 + R;
  } else if (axis=='z') {
    xlo = c1 - R;
    xhi = c1 + R;
    ylo = c2 - R;
    yhi = c2 + R;
    zlo = lo;
    zhi = hi;
  }
}


RegCylinder::~RegCylinder()
{

}

/* ----------------------------------------------------------------------
   inside = 1 if x,y,z is inside or on surface
   inside = 0 if x,y,z is outside and not on surface
------------------------------------------------------------------------- */

int RegCylinder::inside(double x, double y, double z)
{
  //cout << "Check if point (" << x << ", " << y << ", " << z << ") is inside the region" << endl;
  double dSq;
  if (axis=='x') {
    dSq = square(y - c1) + square(z - c2);
    if (x >= lo && x <= hi && dSq <= RSq) {
      return 1;
    }
  } else if (axis=='y') {
    dSq = square(x - c1) + square(z - c2);
    if (y >= lo && y <= hi && dSq <= RSq) {
      return 1;
    }
  } else if (axis=='z') {
    dSq = square(x - c1) + square(y - c2);
    if (z >= lo && z <= hi && dSq <= RSq) {
      return 1;
    }
  }
  return 0;
}

/* ----------------------------------------------------------------------
   return a vector that contains the limits of the box
------------------------------------------------------------------------- */
vector<double> RegCylinder::limits(){
  vector<double> lim;
  lim.push_back(xlo);
  lim.push_back(xhi);
  lim.push_back(ylo);
  lim.push_back(yhi);
  lim.push_back(zlo);
  lim.push_back(zhi);
  return lim;
}
