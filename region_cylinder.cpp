#include <iostream>
#include "domain.h"
#include "region_cylinder.h"
#include "input.h"
#include "var.h"
#include "math_special.h"

using namespace std;
using namespace MathSpecial;

#define BIG 1.0e20

RegCylinder::RegCylinder(MPM *mpm, vector<string> args) : Region(mpm, args)
{
  cout << "Initiate RegCylinder" << endl;

  if (args.size()<8) {
    cout << "Error: region command not enough arguments" << endl;
    exit(1);
  }

  options(&args, args.begin()+8);

  if (args[2].compare("x") == 0) {
    axis = 'x';
  } else if (args[2].compare("y") == 0) {
    axis = 'y';
  } else if (args[2].compare("z") == 0) {
    axis = 'z';
  } else {
    cout << "Error: region cylinder axis not understood, expect x, y, or z, received " << args[2] << endl;
    exit(1);
  }

  c1 = input->parsev(args[3]);
  c2 = input->parsev(args[4]);
  R = input->parsev(args[5]);
  RSq = R*R;

  cout << "axis, c1, c2, R = " << axis << "\t" << c1 << "\t" << c2 << "\t" << R << endl;

  if (args[6].compare("INF") == 0 || args[6].compare("EDGE") == 0) {
    if (domain->regions.size() == 0) {
      cout << "Cannot use region INF or EDGE when box does not exist" << endl;
      exit(1);
    }
    lo = -BIG;
  } else lo = input->parsev(args[6]);

  if (args[7].compare("INF") == 0 || args[7].compare("EDGE") == 0) {
    if (domain->regions.size() == 0) {
      cout << "Cannot use region INF or EDGE when box does not exist" << endl;
      exit(1);
    }
    hi = BIG;
  } else hi = input->parsev(args[7]);

  cout << "lo hi = " << lo << "\t" << hi << endl;

  // error check

  if (lo > hi) {
    cout << "Illegal region cylinder command: low is higher than high" << endl;
    exit(1);
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

  if (domain->boxlo[0] > xlo) domain->boxlo[0] = xlo;
  if (domain->boxlo[1] > ylo) domain->boxlo[1] = ylo;
  if (domain->boxlo[2] > zlo) domain->boxlo[2] = zlo;

  if (domain->boxhi[0] < xhi) domain->boxhi[0] = xhi;
  if (domain->boxhi[1] < yhi) domain->boxhi[1] = yhi;
  if (domain->boxhi[2] < zhi) domain->boxhi[2] = zhi;
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
