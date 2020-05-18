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

#include "region_block.h"
#include "domain.h"
#include "error.h"
#include "input.h"
#include "var.h"
#include <iostream>

using namespace std;

#define BIG 1.0e20

Block_::Block_(MPM *mpm, vector<string> args) : Region(mpm, args)
{
  cout << "Initiate Block_" << endl;

  if (domain->dimension == 3 && args.size()<8) {
    error->all(FLERR, "Error: region command not enough arguments.\n");
  } else if (domain->dimension == 2 && args.size()<6) {
    error->all(FLERR, "Error: region command not enough arguments.\n");
  } else if (domain->dimension == 1 && args.size()<4) {
    error->all(FLERR, "Error: region command not enough arguments.\n");
  }

  if (domain->dimension == 3) options(&args, args.begin()+8);
  else if (domain->dimension == 2) options(&args, args.begin()+6);
  else if (domain->dimension == 1) options(&args, args.begin()+4);

  if (args[2].compare("INF") == 0 || args[2].compare("-INF") == 0 || args[2].compare("EDGE") == 0) {
    if (domain->regions.size() == 0) {
      error->all(FLERR, "Cannot use region INF or EDGE when box does not exist.\n");
    }
    xlo = -BIG;
  } else {
    xlo = input->parsev(args[2]);
    if (domain->boxlo[0] > xlo) domain->boxlo[0] = xlo;
  }

  if (args[3].compare("INF") == 0 || args[3].compare("+INF") == 0 || args[3].compare("EDGE") == 0) {
    if (domain->regions.size() == 0) {
      error->all(FLERR, "Cannot use region INF or EDGE when box does not exist.\n");
    }
    xhi = BIG;
  } else {
    xhi = input->parsev(args[3]);
    if (domain->boxhi[0] < xhi) domain->boxhi[0] = xhi;
  }

  cout << "xlo xhi = " << xlo << "\t" << xhi << endl;

  if (domain->dimension >= 2) {
    if (args[4].compare("-INF") == 0 || args[4].compare("INF") == 0 || args[4].compare("EDGE") == 0) {
      if (domain->regions.size() == 0) {
	error->all(FLERR, "Cannot use region INF or EDGE when box does not exist.\n");
      }
      ylo = -BIG;
    } else {
      ylo = input->parsev(args[4]);
      if (domain->boxlo[1] > ylo) domain->boxlo[1] = ylo;
    }

    if (args[5].compare("INF") == 0 || args[5].compare("+INF") == 0 || args[5].compare("EDGE") == 0) {
      if (domain->regions.size() == 0) {
	error->all(FLERR, "Cannot use region INF or EDGE when box does not exist.\n");
      }
      yhi = BIG;
    } else {
      yhi = input->parsev(args[5]);
      if (domain->boxhi[1] < yhi) domain->boxhi[1] = yhi;
    }

    cout << "ylo yhi = " << ylo << "\t" << yhi << endl;
  } else {
    ylo = yhi = 0;
    zlo = zhi = 0;
  }

  if (domain->dimension == 3) {
    if (args[6].compare("+INF") == 0 || args[6].compare("INF") == 0 || args[6].compare("EDGE") == 0) {
      if (domain->regions.size() == 0) {
	error->all(FLERR, "Cannot use region INF or EDGE when box does not exist.\n");
      }
      zlo = -BIG;
    } else {
      zlo = input->parsev(args[6]);
      if (domain->boxlo[2] > zlo) domain->boxlo[2] = zlo;
    }

    if (args[7].compare("+INF") == 0 || args[7].compare("INF") == 0 || args[7].compare("EDGE") == 0) {
      if (domain->regions.size() == 0) {
	error->all(FLERR, "Cannot use region INF or EDGE when box does not exist.\n");
      }
      zhi = BIG;
    } else {
      zhi = input->parsev(args[7]);
      if (domain->boxhi[2] < zhi) domain->boxhi[2] = zhi;
    }

    cout << "zlo zhi = " << zlo << "\t" << zhi << endl;
  } else {
    zlo = zhi = 0;
  }

  // error check

  if (xlo > xhi || ylo > yhi || zlo > zhi) {
    error->all(FLERR, "Illegal region block command.\n");
  }
}


Block_::~Block_()
{

}

/* ----------------------------------------------------------------------
   inside = 1 if x,y,z is inside or on surface
   inside = 0 if x,y,z is outside and not on surface
------------------------------------------------------------------------- */

int Block_::inside(double x, double y, double z)
{
  //cout << "Check if point (" << x << ", " << y << ", " << z << ") is inside the region" << endl;
  if (x >= xlo && x <= xhi && y >= ylo && y <= yhi && z >= zlo && z <= zhi)
    return 1;
  return 0;
}

/* ----------------------------------------------------------------------
   return a vector that contains the limits of the box
------------------------------------------------------------------------- */
vector<double> Block_::limits(){
  vector<double> lim;
  lim.push_back(xlo);
  lim.push_back(xhi);
  lim.push_back(ylo);
  lim.push_back(yhi);
  lim.push_back(zlo);
  lim.push_back(zhi);
  return lim;
}
