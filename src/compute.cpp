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

#include "compute.h"
#include "group.h"
#include "input.h"
#include "var.h"
#include <iostream>

using namespace std;


Compute::Compute(MPM *mpm, vector<string> args) :
  Pointers(mpm)
{
  cout << "Creating new compute with ID: " << args[0] << endl;
  id = args[0];

  igroup = group->find(args[2]);
  if (igroup == -1) {
    cout << "Could not find group ID " << args[2] << endl;
  }
  groupbit = group->bitmask[igroup];
}
