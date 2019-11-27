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

#include <iostream>
#include <algorithm>
#include "dump.h"
#include "error.h"

using namespace std;


Dump::Dump(MPM *mpm, vector<string> args) :
  Pointers(mpm)
{
  cout << "Creating new dump with ID: " << args[0] << endl;
  id = args[0];

  if (args[1].compare("all")!=0) {
    error->all(FLERR, "Error: groups are not yet supported for dumps. Should use all.\n");
  }

  style = args[2];
  
  filename = args[4];
}

Dump::~Dump()
{
}
