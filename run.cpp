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

/* ---------------------------------------------------------------------- */

Run::Run(MPM *mpm) : Pointers(mpm) {}

/* ---------------------------------------------------------------------- */

Var Run::command(vector<string> args)
{
  // cout << "In Run::command()" << endl;

  if (args.size() < 1) {
    cout << "Illegal run command" << endl;
    exit(1);
  }

  mpm->init();

  update->scheme->setup();

  // Check that a method is available:
  if (update->method == NULL) {
    cout << "Error: no method was defined!" << endl;
    exit(1);
  }

  int nsteps = (int) input->parsev(args[0]);
  update->nsteps = nsteps;//INT_MAX;
  update->maxtime = -1;
  update->firststep = update->ntimestep;
  update->laststep = update->firststep + nsteps;
  update->scheme->run(Var("timestep<"+to_string(update->laststep), update->ntimestep<update->laststep));

  return Var(0);
}
