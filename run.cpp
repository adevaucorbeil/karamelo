#include <iostream>
#include <vector>
#include "run.h"
#include "domain.h"
#include "output.h"
#include "input.h"
#include "update.h"
#include "scheme.h"
#include "var.h"

/* ---------------------------------------------------------------------- */

Run::Run(MPM *mpm) : Pointers(mpm) {}

/* ---------------------------------------------------------------------- */

void Run::command(vector<string> args)
{
  cout << "In Run::command()" << endl;

  if (args.size() < 1) {
    cout << "Illegal run command" << endl;
    exit(1);
  }

  update->scheme->setup();

  int nsteps = (int) input->parsev(args[0]);
  update->nsteps = nsteps;
  update->firststep = update->ntimestep + 1;
  update->laststep = update->firststep + nsteps;
  update->scheme->run(nsteps);
}
