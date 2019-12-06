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

#include "fix_max_plastic_strain.h"
#include "domain.h"
#include "error.h"
#include "group.h"
#include "input.h"
#include "math_special.h"
#include "output.h"
#include "update.h"
#include <Eigen/Eigen>
#include <iostream>
#include <vector>

using namespace std;
using namespace FixConst;
using namespace MathSpecial;
using namespace Eigen;

FixMaxPlasticStrain::FixMaxPlasticStrain(MPM *mpm, vector<string> args)
    : Fix(mpm, args) {
  if (args.size() < 3) {
    error->all(FLERR, "Error: too few arguments for fix_max_plastic_strain: "
                      "requires at least 3 arguments. " +
                          to_string(args.size()) + " received\n");
  }

  if (group->pon[igroup].compare("particles") != 0 &&
      group->pon[igroup].compare("all") != 0) {
    error->all(
        FLERR,
        "fix_max_plastic_strain needs to be given a group of particles " +
            group->pon[igroup] + ", " + args[2] + " is a group of " +
            group->pon[igroup] + ".\n");
  }

  cout << "Creating new fix FixMaxPlasticStrain with ID: " << args[0] << endl;
  id = args[0];

  (*input->vars)[id + "_1s"] = Var(id + "_1s", 0);
  (*input->vars)[id + "_2s"] = Var(id + "_2s", 0);
}

FixMaxPlasticStrain::~FixMaxPlasticStrain() {}

void FixMaxPlasticStrain::init() {}

void FixMaxPlasticStrain::setup() {}

void FixMaxPlasticStrain::setmask() {
  mask = 0;
  mask |= FINAL_INTEGRATE;
}

void FixMaxPlasticStrain::final_integrate() {
  if (update->ntimestep != output->next && update->ntimestep != update->nsteps)
    return;
  // cout << "In FixMaxPlasticStrain::post_particles_to_grid()\n";

  // Go through all the nodes in the group and set b to the right value:
  double ux, uy, uz;

  int solid = group->solid[igroup];

  Solid *s;

  double Es(0.), Tmax(0.);

  if (solid == -1) {
    for (int isolid = 0; isolid < domain->solids.size(); isolid++) {
      s = domain->solids[isolid];

      for (int in = 0; in < s->np_local; in++) {
        if (s->mask[in] & groupbit) {
          Es = max(Es, s->eff_plastic_strain[in]);
          Tmax = max(Tmax, s->T[in]);
        }
      }
    }

    (*input->vars)[id + "_1s"] = Var(id + "_1s", Es);
    (*input->vars)[id + "_2s"] = Var(id + "_2s", Tmax);
  } else {
    s = domain->solids[solid];

    for (int in = 0; in < s->np_local; in++) {
      if (s->mask[in] & groupbit) {
	Es = max(Es, s->eff_plastic_strain[in]);
	Tmax = max(Tmax, s->T[in]);
      }
    }

    (*input->vars)[id + "_1s"] = Var(id + "_1s", Es);
    (*input->vars)[id + "_2s"] = Var(id + "_2s", Tmax);
  }
}
