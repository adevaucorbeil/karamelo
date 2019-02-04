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

  xlo = input->parse(args[2]);
  xhi = input->parse(args[3]);
  cout << "xlo xhi = " << xlo << "\t" << xhi << endl;
  ylo = input->parse(args[4]);
  yhi = input->parse(args[5]);
  cout << "ylo yhi = " << ylo << "\t" << yhi << endl;
  zlo = input->parse(args[6]);
  zhi = input->parse(args[7]);
  cout << "zlo zhi = " << zlo << "\t" << zhi << endl;
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
  cout << "Check if point (" << x << ", " << y << ", " << z << ") is inside the region" << endl;
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
