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

#include <iostream>
#include <vector>
#include <climits>
#include <run_time.h>
#include <domain.h>
#include <output.h>
#include <input.h>
#include <update.h>
#include <scheme.h>
#include <var.h>
#include <error.h>

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

  float maxtime = (float) input->parsev(args[0]) + update->atime;
  update->nsteps = INT_MAX;
  update->maxtime = maxtime;
  update->firststep = update->ntimestep;
  update->laststep = INT_MAX;
  update->scheme->run(Var("time<"+to_string(maxtime), update->atime < maxtime ));

  return Var(0);
}
