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
#include "compute.h"

using namespace std;


Compute::Compute(MPM *mpm, vector<string> args) :
  Pointers(mpm)
{
  cout << "Creating new compute with ID: " << args[0] << endl;
  id = args[0];
}
