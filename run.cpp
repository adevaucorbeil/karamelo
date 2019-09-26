#include <iostream>
#include <vector>
#include "run.h"
#include "domain.h"
#include "output.h"
#include "input.h"
#include "update.h"
#include "scheme.h"
#include "var.h"
#include <string>
#include "error.h"

/* ---------------------------------------------------------------------- */

Run::Run(MPM *mpm) : Pointers(mpm) {}

/* ---------------------------------------------------------------------- */

Var Run::command(vector<string> args)
{
  // cout << "In Run::command()" << endl;

  if (args.size() < 1) {
    error->all(FLERR, "Illegal run command.\n");
  }

  mpm->init();

  update->scheme->setup();

  // Check that a method is available:
  if (update->method == NULL) {
    error->all(FLERR, "Error: no method was defined!\n");
  }

  int nsteps = (int) input->parsev(args[0]);
  update->nsteps = INT_MAX;
  update->maxtime = -1;
  update->firststep = update->ntimestep;
  update->laststep = update->firststep + nsteps;
  update->scheme->run(Var("timestep<"+to_string(update->laststep), update->ntimestep<update->laststep));

  return Var(0);
}
