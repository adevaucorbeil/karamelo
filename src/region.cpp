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

#include "region.h"
#include "error.h"
#include "universe.h"
#include <iostream>

using namespace std;

Region::Region(MPM *mpm, vector<string> args) : Pointers(mpm)
{
  if (universe->me == 0)
    cout << "Creating new region with ID: " << args[0] << endl;
  id       = args[0];
  style    = args[1];
  interior = 1; // Interior by default
}
Region::~Region() {}

void Region::init() {}

void Region::options(vector<string> *args, vector<string>::iterator it)
{

  //cout << "In region::options()" << endl;
  if (args->end() < it)
  {
    // error->all(FLERR, "Error: not enough arguments.\n");
    return;
  }
  if (args->end() > it)
  {
    for (it; it != args->end(); ++it)
    {
      if ((*it).compare("exterior") == 0)
	{
	  if (universe->me == 0)
	    cout << "\nRecognized exterior argument\n";
	  interior = 0;
	}
      else
	{
	  if (universe->me == 0) {
	    cout << "Ignoring optional arguments: ";
	    cout << *it << "\t";
	  }
	}
      if (universe->me == 0)
	cout << endl;
    }
  }
}

int Region::match(double x, double y, double z)
{
  if (interior)
    return inside(x, y, z);
  else
  {
    return !(inside(x, y, z));
  }
  // return !(inside(x,y,z) ^ interior);
}
