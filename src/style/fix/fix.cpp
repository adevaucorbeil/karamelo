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

Fix::Fix(MPM *mpm, const vector<string> &args, int mask):
  Pointers(mpm),
  id(args.at(0)),
  style(args.at(1)),
  mask(mask),
  igroup(group->find(args.at(2)))
{
  if (!universe->me)
  {
    cout << "Creating new fix with ID: " << args[0] << endl;

    if (igroup == -1)
      cout << "Could not find group ID " << args[2] << endl;
  }

  groupbit = group->bitmask[igroup];

  (*input->vars)[id + "_x"] = Var(id + "_x", 0);
  (*input->vars)[id + "_y"] = Var(id + "_y", 0);
  (*input->vars)[id + "_z"] = Var(id + "_z", 0);
  (*input->vars)[id + "_s"] = Var(id + "_s", 0);
}
