#include <iostream>
#include <vector>
#include "run_time.h"
#include "domain.h"
#include "output.h"
#include "input.h"
#include "update.h"
#include "scheme.h"
#include "var.h"

/* ---------------------------------------------------------------------- */

RunTime::RunTime(MPM *mpm) : Pointers(mpm) {}

/* ---------------------------------------------------------------------- */

Var RunTime::command(vector<string> args)
{
  // cout << "In RunTime::command()" << endl;

  if (args.size() < 1) {
    cout << "Illegal run command" << endl;
    exit(1);
  }

  mpm->init();

  update->scheme->setup();

  double maxtime = (double) input->parsev(args[0]);
  update->nsteps = INT_MAX;
  update->maxtime = maxtime;
  update->firststep = update->ntimestep + 1;
  update->laststep = INT_MAX;
  update->scheme->run(update->nsteps);

  return Var(0);
}
