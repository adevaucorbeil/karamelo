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
#include <vector>
#include "write_restart.h"
#include "domain.h"
#include "output.h"
#include "input.h"
#include "update.h"
#include "scheme.h"
#include "var.h"
#include "error.h"

/* ---------------------------------------------------------------------- */

WriteRestart::WriteRestart(MPM *mpm) : Pointers(mpm) {}

Var WriteRestart::command(vector<string> args)
{
  cout << "In WriteRestart::command()" << endl;

  if (args.size() < 1) {
    error->all(FLERR, "Illegal write command.\n");
  }

  return Var(0);
}


void WriteRestart::write(){
  cout << "In WriteRestart::write()" << endl;
}
