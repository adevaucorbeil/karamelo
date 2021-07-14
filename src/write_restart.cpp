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
#include "domain.h"
#include "error.h"
#include "group.h"
#include "method.h"
#include "modify.h"
#include "universe.h"
#include "update.h"
#include "var.h"
#include "version.h"
#include <iostream>
#include <vector>

enum { VERSION, DIMENSION, NPROCS };
/* ---------------------------------------------------------------------- */

WriteRestart::WriteRestart(MPM *mpm) : Pointers(mpm) {
  of = NULL;
}

WriteRestart::~WriteRestart() {
  delete of;
}

Var WriteRestart::command(vector<string> args) {
  // cout << "In WriteRestart::command()" << endl;

  if (args.size() < 1) {
    error->all(FLERR, "Illegal write command.\n");
  }

  filename = args[0];
  pos_asterisk = filename.find('*'); // Check if there is a * in the name, and record its position.
  return Var(0);
}


void WriteRestart::write() {

  string frestart;
  if (pos_asterisk != string::npos) {
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
    frestart = filename + "proc-" + to_string(universe->me) + ".";

  if (universe->me == 0)
    cout << "write " << frestart << endl;
  of = new ofstream(frestart, ios_base::out | ios_base::binary);

  if (of->is_open()) {

    // proc 0 writes out header:
    if (universe->me == 0)
      header();

    // Everyone writes the method, scheme, timestep, dt:
    update->write_restart(of);
    domain->write_restart(of);
    group->write_restart(of);
    modify->write_restart(of);


    of->close();
  } else {
    error->all(FLERR, "Error: cannot write in file: " + frestart + ".\n");
  }
}

/*! Writes out problem description in restart file.
 */
void WriteRestart::header() {
  // Karamelo version:
  write_string(VERSION, Version::GIT_SHA1);
  write_variable(DIMENSION, domain->dimension);
  write_variable(NPROCS, universe->nprocs);

  // -1 flag signals end of header
  int flag = -1;
  of->write(reinterpret_cast<const char *>(&flag), sizeof(int));
}

/*!  write a flag and a variable (which is not a string) into restart file.
 */
template <typename T> void WriteRestart::write_variable(int flag, T value) {
  of->write(reinterpret_cast<const char *>(&flag), sizeof(int));
  of->write(reinterpret_cast<const char *>(&value), sizeof(value));
}


/*!  write a flag and a string into restart file.
 */
void WriteRestart::write_string(int flag, string value) {
  size_t N = value.size();
  of->write(reinterpret_cast<const char *>(&flag), sizeof(int));
  of->write(reinterpret_cast<const char *>(&N), sizeof(size_t));
  of->write(reinterpret_cast<const char *>(value.c_str()), N);
}

// /*!  write a flag and a variable into restart file.
//  */
// void WriteRestart::write_int(int flag, int value) {
//   of->write(reinterpret_cast<const char *>(&flag), sizeof(int));
//   of->write(reinterpret_cast<const char *>(&value), sizeof(int));
// }
