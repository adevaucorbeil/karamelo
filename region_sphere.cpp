#include <iostream>
#include "domain.h"
#include "region_sphere.h"
#include "input.h"
#include "var.h"
#include "math_special.h"

using namespace std;
using namespace MathSpecial;

#define BIG 1.0e20

RegSphere::RegSphere(MPM *mpm, vector<string> args) : Region(mpm, args)
{
  cout << "Initiate RegSphere" << endl;

  c1 = c2 = c3 = 0;
  R = 0;

  if (args.size() < Nargs[domain->dimension-1]) {
    cout << "Error: region command not enough arguments" << endl;
    cout << usage[domain->dimension-1];
    exit(1);
  }
  
  if (args.size() > Nargs[domain->dimension-1]) {
    cout << "Error: region command too many arguments" << endl;
    cout << usage[domain->dimension-1];
    exit(1);
  }

  int iargs=2;
  c1 = input->parsev(args[iargs++]);
  if (domain->dimension >= 2) c2 = input->parsev(args[iargs++]);
  if (domain->dimension == 3) c3 = input->parsev(args[iargs++]);
  R = input->parsev(args[iargs++]);

  RSq = R*R;
  cout << "c1, c2, c3, R = " << c1 << "\t" << c2 << "\t" << c3 << "\t" << R << endl;

  xlo = c1 - R;
  xhi = c1 + R;
  ylo = c2 - R;
  yhi = c2 + R;
  zlo = c3 - R;
  zhi = c3 + R;

  if (domain->boxlo[0] > xlo) domain->boxlo[0] = xlo;
  if (domain->boxlo[1] > ylo) domain->boxlo[1] = ylo;
  if (domain->dimension == 3) if (domain->boxlo[2] > zlo) domain->boxlo[2] = zlo;

  if (domain->boxhi[0] < xhi) domain->boxhi[0] = xhi;
  if (domain->boxhi[1] < yhi) domain->boxhi[1] = yhi;
  if (domain->dimension == 3) if (domain->boxhi[2] < zhi) domain->boxhi[2] = zhi;
}


RegSphere::~RegSphere()
{

}

/* ----------------------------------------------------------------------
   inside = 1 if x,y,z is inside or on surface
   inside = 0 if x,y,z is outside and not on surface
------------------------------------------------------------------------- */

int RegSphere::inside(double x, double y, double z)
{
  //cout << "Check if point (" << x << ", " << y << ", " << z << ") is inside the region" << endl;
  double dSq;
  dSq = square(x - c1) + square(y - c2) + square(z - c3);
  if (dSq <= RSq) {
    return 1;
  }
  return 0;
}

/* ----------------------------------------------------------------------
   return a vector that contains the limits of the box
------------------------------------------------------------------------- */
vector<double> RegSphere::limits(){
  vector<double> lim;
  lim.push_back(xlo);
  lim.push_back(xhi);
  lim.push_back(ylo);
  lim.push_back(yhi);
  lim.push_back(zlo);
  lim.push_back(zhi);
  return lim;
}
