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

#include "read_restart.h"
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

ReadRestart::ReadRestart(MPM *mpm) : Pointers(mpm) {
  ifr = NULL;
}

ReadRestart::~ReadRestart() {
  delete ifr;
}

Var ReadRestart::command(vector<string> args) {
  cout << "In ReadRestart::command()" << endl;

  if (args.size() < 1) {
    error->all(FLERR, "Illegal read command.\n");
  }

  filename = args[0];
  pos_asterisk = filename.find('%'); // Check if there is a * in the name, and record its position.

  string frestart;
  if (pos_asterisk != string::npos) {
    // Replace the asterisk by proc-N.ntimestep:
    frestart = filename.substr(0, pos_asterisk);
    if (universe->nprocs > 1) {
      frestart += "proc-" + to_string(universe->me) + ".";
    }
    if (filename.size() - pos_asterisk - 1 > 0)
      frestart +=
          filename.substr(pos_asterisk + 1, filename.size() - pos_asterisk - 1);
  } else
    frestart = filename;

  cout << "read " << frestart << endl;
  ifr = new ifstream(frestart, ios_base::in | ios_base::binary);

  if (ifr->is_open()) {

    // proc 0 reads out header:
    if (universe->me == 0)
      header();

    MPI_Bcast(&domain->dimension, 1, MPI_INT, 0, universe->uworld);

    // Everyone reads the method, scheme, timestep, dt:
    update->read_restart(ifr);
    domain->read_restart(ifr);
    group->read_restart(ifr);
    modify->read_restart(ifr);


    ifr->close();
  } else {
    error->all(FLERR, "Error: cannot read in file: " + frestart + ".\n");
  }
  return Var(0);
}

/*! Reads out problem description in restart file.
 */
void ReadRestart::header() {
  // Karamelo version:
  int flag = read_int();

  while (flag >= 0) {
    if (flag == VERSION) {
      string version = read_string();
      cout << "version = " << version << endl;
    } else if (flag == DIMENSION) {
      domain->dimension = read_int();
      cout << "dimension = " << domain->dimension << endl;
    } else if (flag == NPROCS) {
      int nprocs = read_int();
      cout << "nprocs = " << nprocs << endl;
      if (nprocs != universe->nprocs) {
	error->one(FLERR, "Restart file written for " + to_string(nprocs) + " CPUs.\n");
      }
    }
    flag = read_int();
  }
}

/*!  read a flag and a variable into restart file.
 */
int ReadRestart::read_int() {

  int i = -1;
  ifr->read(reinterpret_cast<char *>(&i), sizeof(int));
  return i;
}


/*!  read a flag and a variable into restart file.
 */
string ReadRestart::read_string() {

  size_t N = 0;
  string s = "";

  ifr->read(reinterpret_cast<char *>(&N), sizeof(size_t));
  s.resize(N);

  ifr->read(reinterpret_cast<char *>(&s[0]), N);
  return s;
}
