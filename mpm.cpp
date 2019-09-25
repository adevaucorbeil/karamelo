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

MPM::MPM(int narg, char **arg, MPI_Comm communicator)
{
  memory = new Memory(this);
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
  while (iarg < narg) {
    if (strcmp(arg[iarg],"-in") == 0 ||
	strcmp(arg[iarg],"-i") == 0) {
      if (iarg+2 > narg) {
	printf("Invalid command-line argument\n");
	exit(1);
      }
      inflag = iarg + 1;
      iarg += 2;
    }
    iarg++;
  }

  if (universe->me == 0) {
    if (inflag != 0) infile.open(arg[inflag], ios_base::in); // open in read only
    //logfile.open("log.mpm", ios_base::out); // open in write only
    wlogfile = new ofstream("log.mpm", ios_base::out);

    if (!wlogfile->is_open()){
      printf("Cannot open file log.mpm\n");
      exit(1);
    }
    if (!infile.is_open()) {
      printf("Cannot open input script %s\n",arg[inflag]);
      exit(1);
    }
  }
}

MPM::~MPM()
{
  delete memory;
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

  if (wlogfile) wlogfile = NULL;

  delete variables;
}

void MPM::init() {
  modify->init();
}
