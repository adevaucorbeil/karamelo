#include <iostream>
#include "domain.h"
#include "region_block.h"
#include "input.h"

using namespace std;

#define BIG 1.0e20

RegBlock::RegBlock(MPM *mpm, vector<string> args) : Region(mpm, args)
{
  cout << "Initiate RegBlock" << endl;

  if (args.size()<8) {
    cout << "Error: region command not enough arguments" << endl;
    exit(1);
  }

  options(&args, args.begin()+8);

  if (args[2].compare("INF") == 0 || args[2].compare("EDGE") == 0) {
    if (domain->regions.size() == 0) {
      cout << "Cannot use region INF or EDGE when box does not exist" << endl;
      exit(1);
    }
    if (args[2].compare("INF") == 0) xlo = -BIG;
  } else xlo = input->parse(args[2]);

  if (args[3].compare("INF") == 0 || args[3].compare("EDGE") == 0) {
    if (domain->regions.size() == 0) {
      cout << "Cannot use region INF or EDGE when box does not exist" << endl;
      exit(1);
    }
    if (args[3].compare("INF") == 0) xhi = BIG;
  } else xhi = input->parse(args[3]);

  cout << "xlo xhi = " << xlo << "\t" << xhi << endl;

  if (args[4].compare("INF") == 0 || args[4].compare("EDGE") == 0) {
    if (domain->regions.size() == 0) {
      cout << "Cannot use region INF or EDGE when box does not exist" << endl;
      exit(1);
    }
    if (args[4].compare("INF") == 0) ylo = -BIG;
  } else ylo = input->parse(args[4]);

  if (args[5].compare("INF") == 0 || args[5].compare("EDGE") == 0) {
    if (domain->regions.size() == 0) {
      cout << "Cannot use region INF or EDGE when box does not exist" << endl;
      exit(1);
    }
    if (args[5].compare("INF") == 0) yhi = BIG;
  } else yhi = input->parse(args[5]);

  cout << "ylo yhi = " << ylo << "\t" << yhi << endl;

  if (args[6].compare("INF") == 0 || args[6].compare("EDGE") == 0) {
    if (domain->regions.size() == 0) {
      cout << "Cannot use region INF or EDGE when box does not exist" << endl;
      exit(1);
    }
    if (args[6].compare("INF") == 0) zlo = -BIG;
  } else zlo = input->parse(args[6]);

  if (args[7].compare("INF") == 0 || args[7].compare("EDGE") == 0) {
    if (domain->regions.size() == 0) {
      cout << "Cannot use region INF or EDGE when box does not exist" << endl;
      exit(1);
    }
    if (args[7].compare("INF") == 0) zhi = BIG;
  } else yhi = input->parse(args[7]);

  cout << "zlo zhi = " << zlo << "\t" << zhi << endl;

  // error check

  if (xlo > xhi || ylo > yhi || zlo > zhi) {
    cout << "Illegal region block command" << endl;
    exit(1);
  }

  if (domain->boxlo[0] > xlo) domain->boxlo[0] = xlo;
  if (domain->boxlo[1] > ylo) domain->boxlo[1] = ylo;
  if (domain->boxlo[2] > zlo) domain->boxlo[2] = zlo;

  if (domain->boxhi[0] < xhi) domain->boxhi[0] = xhi;
  if (domain->boxhi[1] < yhi) domain->boxhi[1] = yhi;
  if (domain->boxhi[2] < zhi) domain->boxhi[2] = zhi;
}


RegBlock::~RegBlock()
{

}

/* ----------------------------------------------------------------------
   inside = 1 if x,y,z is inside or on surface
   inside = 0 if x,y,z is outside and not on surface
------------------------------------------------------------------------- */

int RegBlock::inside(double x, double y, double z)
{
  //cout << "Check if point (" << x << ", " << y << ", " << z << ") is inside the region" << endl;
  if (x >= xlo && x <= xhi && y >= ylo && y <= yhi && z >= zlo && z <= zhi)
    return 1;
  return 0;
}

/* ----------------------------------------------------------------------
   return a vector that contains the limits of the box
------------------------------------------------------------------------- */
vector<double> RegBlock::limits(){
  vector<double> lim;
  lim.push_back(xlo);
  lim.push_back(xhi);
  lim.push_back(ylo);
  lim.push_back(yhi);
  lim.push_back(zlo);
  lim.push_back(zhi);
  return lim;
}
