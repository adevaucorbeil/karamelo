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
#include <group.h>
#include <var.h>
#include <internal_force.h>
#include <error.h>

using namespace std;

IntForce::IntForce(MPM *mpm) : Pointers(mpm) {}

Var IntForce::command(vector<string> args) {
  // cout << "In IntForce::command()" << endl;

  if (args.size() < 2) {
    error->all(FLERR,"Illegal run command.\n");
  }

  igroup = group->find(args[0]);

  if (igroup == -1) {
    error->all(FLERR, "Error: could not find group named: " + args[0] + ".\n");
  }

  if (args[1] == "x") dir = 0;
  else if (args[1] == "y") dir = 1;
  else if (args[1] == "z") dir = 2;
  else {
    error->all(FLERR, "Error: directions should be either x,y or z: " + args[1] + " not understood.\n");
  }

  float iforce = group->internal_force(igroup, dir);
  string eq = "internal_force(" + args[0] + "," + args[1] + ")";
  return Var(eq, iforce);
}
