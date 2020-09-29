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

#include "write_restart.h"
#include "error.h"
#include "universe.h"
#include "update.h"
#include "var.h"
#include <iostream>
#include <vector>

/* ---------------------------------------------------------------------- */

WriteRestart::WriteRestart(MPM *mpm) : Pointers(mpm) {}

Var WriteRestart::command(vector<string> args) {
  cout << "In WriteRestart::command()" << endl;

  if (args.size() < 1) {
    error->all(FLERR, "Illegal write command.\n");
  }

  filename = args[0];
  pos_asterisk = filename.find('*'); // Check if there is a * in the name, and record its position.

  return Var(0);
}

void WriteRestart::write() {

  string frestart;
  if (pos_asterisk >= 0) {
    // Replace the asterisk by proc-N.ntimestep:
    frestart = filename.substr(0, pos_asterisk);
    if (universe->nprocs > 1) {
      frestart += "proc-" + to_string(universe->me) + ".";
    }
    frestart += to_string(update->ntimestep);
    if (filename.size() - pos_asterisk - 1 > 0)
      frestart +=
          filename.substr(pos_asterisk + 1, filename.size() - pos_asterisk - 1);
  } else
    frestart = filename;


  cout << "write " << frestart << endl;
}
