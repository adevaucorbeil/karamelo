#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mpm.h"
#include "domain.h"
#include "input.h"

MPM::MPM(int narg, char **arg)
{
  input = new Input(this, narg, arg);

  domain = new Domain(this);
  //settings = new Settings(this);

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
  }

  if (inflag != 0) infile.open(arg[inflag], ios_base::in); // open in read only

  if (!infile.is_open()) {
    printf("Cannot open input script %s\n",arg[inflag]);
    exit(1); 
  }
}

MPM::~MPM()
{
  delete input;

  if (infile.is_open()) infile.close();
  if (logfile.is_open()) logfile.close();
}
