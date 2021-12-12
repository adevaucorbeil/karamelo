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

#include <fix.h>
#include <group.h>
#include <input.h>
#include <universe.h>
#include <var.h>
#include <iostream>

using namespace std;


Fix::Fix(MPM *mpm, vector<string> args_) :
  Pointers(mpm)
{
  mask = 0;
  args = args_;
  if (universe->me == 0) {
    cout << "Creating new fix with ID: " << args[0] << endl;
  }
  id = args[0];

  style = args[1];

  igroup = group->find(args[2]);
  if (igroup == -1 && universe->me == 0) {
    cout << "Could not find group ID " << args[2] << endl;
  }
  groupbit = group->bitmask[igroup];
  (*input->vars)[id+"_x"]=Var(id+"_x", 0);
  (*input->vars)[id+"_y"]=Var(id+"_y", 0);
  (*input->vars)[id+"_z"]=Var(id+"_z", 0);
  (*input->vars)[id+"_s"]=Var(id+"_s", 0);

}
