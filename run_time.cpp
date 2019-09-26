#include <iostream>
#include <vector>
#include "run_time.h"
#include "domain.h"
#include "output.h"
#include "input.h"
#include "update.h"
#include "scheme.h"
#include "var.h"
#include "error.h"

/* ---------------------------------------------------------------------- */

RunTime::RunTime(MPM *mpm) : Pointers(mpm) {}

/* ---------------------------------------------------------------------- */

Var RunTime::command(vector<string> args)
{
  // cout << "In RunTime::command()" << endl;

  if (args.size() < 1) {
    error->all(FLERR, "Illegal run command.\n");
  }

  mpm->init();

  update->scheme->setup();

  double maxtime = (double) input->parsev(args[0]) + update->atime;
  update->nsteps = INT_MAX;
  update->maxtime = maxtime;
  update->firststep = update->ntimestep;
  update->laststep = INT_MAX;
  update->scheme->run(Var("time<"+to_string(maxtime), update->atime < maxtime ));

  return Var(0);
}
