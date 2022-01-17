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

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <string.h>
#include "mpm.h"
#include "domain.h"
#include "material.h"
#include "input.h"
#include "output.h"
#include "update.h"
#include "modify.h"
#include "memory.h"
#include "group.h"
#include "universe.h"
#include "error.h"
#include "version.h"

MPM::MPM(int narg, char **arg, MPI_Comm communicator)
{
  memory = new Memory(this);
  error = new Error(this);
  universe = new Universe(this, communicator);
  input = new Input(this, narg, arg);
  output = new Output(this);
  update = new Update(this);

  domain = new Domain(this);
  material = new Material(this);
  modify = new Modify(this);
  group = new Group(this);

  initclock = MPI_Wtime();

  // parse input switches

  int inflag = 0;

  // Process command line arguments:
  int iarg = 1;

  if (narg < 2)
    {
      if (universe->me == 0) help();
      error->done(0);
    }

  while (iarg < narg) {
    if (strcmp(arg[iarg],"-in") == 0 ||
	strcmp(arg[iarg],"-i") == 0) {
      if (iarg+2 > narg) {
	error->all(FLERR,"Invalid command-line argument\n");
      }
      inflag = iarg + 1;
      iarg += 2;
    }
    iarg++;
  }

  wlogfile = nullptr;

  if (universe->me == 0) {
    cout << "Karamelo -- Parallel Material Point Methods Simulator Build SHA1:"
         << Version::GIT_SHA1 << endl;
    cout << "Running on " << universe->nprocs << " procs\n";

    if (inflag != 0)
      {
	infile.open(arg[inflag], ios_base::in); // open in read only
      }
    //logfile.open("log.mpm", ios_base::out); // open in write only
    wlogfile = new ofstream("log.mpm", ios_base::out);

    if (!wlogfile->is_open())
      {
	error->one(FLERR,"Cannot open file log.mpm\n");
      }
    if (!infile.is_open())
      {
    	string problem(arg[inflag]);
    	error->one(FLERR,"Cannot open input script " + problem + ".\n");
      }
  }
}

MPM::~MPM()
{
  delete memory;
  delete error;
  delete input;
  delete output;
  delete update;

  delete domain;
  delete material;
  delete modify;
  delete group;

  double totalclock = MPI_Wtime() - initclock;

  int seconds = fmod(totalclock,60.0);
  totalclock  = (totalclock - seconds) / 60.0;
  int minutes = fmod(totalclock,60.0);
  int hours = (totalclock - minutes) / 60.0;
  cout << "Total wall time: " << hours << ":" << minutes << ":" << seconds << endl;

  if (universe->me == 0) {
    if (infile.is_open()) infile.close();
    if (wlogfile->is_open()) wlogfile->close();
  }

  if (wlogfile) wlogfile = nullptr;
}

void MPM::init() {
  modify->init();
}

void MPM::help() {
  // general help message about command line and flags
  cout << "Karamelo -- Parallel Material Point Methods Simulator Build SHA1:"
       << Version::GIT_SHA1 << endl
       << "Usage: karamelo -i input_file\n";
}
