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
#include "group.h"
#include "var.h"
#include "centre_of_mass.h"
#include "error.h"

using namespace std;

CentreOfMass::CentreOfMass(MPM *mpm) : Pointers(mpm) {}

Var CentreOfMass::command(vector<string> args) {
  // cout << "In CentreOfMass::command()" << endl;

  if (args.size() < 2) error->all(FLERR, "Illegal run command");

  igroup = group->find(args[0]);

  if (igroup == -1) {
    error->all(FLERR, "Error: could not find group named: " + args[0] + "\n");
  }

  if (args[1].compare("x") == 0) dir = 0;
  else if (args[1].compare("y") == 0) dir = 1;
  else if (args[1].compare("z") == 0) dir = 2;
  else error->all(FLERR, "Error: directions should be either x,y or z: " + args[1] + " not understood.\n");

  double com = group->xcm(igroup, dir);
  string eq = "xcm(" + args[0] + "," + args[1] + ")";
  return Var(eq, com);
}
