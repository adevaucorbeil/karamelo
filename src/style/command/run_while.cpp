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
#include <run_while.h>
#include <domain.h>
#include <output.h>
#include <input.h>
#include <update.h>
#include <scheme.h>
#include <var.h>
#include <error.h>

/* ---------------------------------------------------------------------- */

RunWhile::RunWhile(MPM *mpm) : Pointers(mpm) {}

/* ---------------------------------------------------------------------- */

Var RunWhile::command(vector<string> args)
{
  // cout << "In RunWhile::command()" << endl;

  if (args.size() < 1) {
    error->all(FLERR, "Illegal run command.\n");
  }

  mpm->init();

  update->scheme->setup();

  Var condition = input->parsev(args[0]);
  update->nsteps = INT_MAX;
  update->maxtime = -1;
  update->firststep = update->ntimestep;
  update->laststep = INT_MAX;
  update->scheme->run(condition);

  return Var(0);
}
