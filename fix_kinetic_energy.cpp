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
#include <string>
#include <Eigen/Eigen>
#include "fix_kinetic_energy.h"
#include "input.h"
#include "group.h"
#include "domain.h"
#include "input.h"
#include "update.h"
#include "output.h"
#include "math_special.h"
#include "universe.h"
#include "error.h"

using namespace std;
using namespace FixConst;
using namespace MathSpecial;
using namespace Eigen;

FixKineticEnergy::FixKineticEnergy(MPM *mpm, vector<string> args) : Fix(mpm, args)
{
  if (args.size() < 3) {
    error->all(FLERR,"Error: too few arguments for fix_kinetic_energy: requires at least 3 arguments. " + to_string(args.size()) + " received.\n");
  }

  if (group->pon[igroup].compare("particles") !=0 && group->pon[igroup].compare("all") !=0) {
    error->all(FLERR, "fix_kinetic_energy needs to be given a group of particles" + group->pon[igroup] + ", " + args[2] + " is a group of " + group->pon[igroup] + ".\n");
  }

  cout << "Creating new fix FixKineticEnergy with ID: " << args[0] << endl;
  id = args[0];
}

FixKineticEnergy::~FixKineticEnergy()
{
}

void FixKineticEnergy::init()
{
}

void FixKineticEnergy::setup()
{
}

void FixKineticEnergy::setmask() {
  mask = 0;
  mask |= FINAL_INTEGRATE;
}


void FixKineticEnergy::final_integrate() {
  if (update->ntimestep != output->next && update->ntimestep != update->nsteps) return;
  // cout << "In FixKineticEnergy::post_particles_to_grid()\n";

  // Go through all the nodes in the group and set b to the right value:
  double ux, uy, uz;
    
  int solid = group->solid[igroup];

  Solid *s;

  double Ek, Ek_reduced;

  Ek = 0;
  Ek_reduced = 0;

  if (solid == -1) {
    for (int isolid = 0; isolid < domain->solids.size(); isolid++) {
      s = domain->solids[isolid];

      for (int in = 0; in < s->np_local; in++) {
	if (s->mask[in] & groupbit) {
	  Ek += 0.5*s->mass[in]*square(s->v[in].norm());
	}
      }
    }
  } else {
    s = domain->solids[solid];

    for (int in = 0; in < s->np_local; in++) {
      if (s->mask[in] & groupbit) {
	Ek += 0.5*s->mass[in]*square(s->v[in].norm());
      }
    }
  }

  // Reduce Ek:
  MPI_Allreduce(&Ek,&Ek_reduced,1,MPI_DOUBLE,MPI_SUM,universe->uworld);
  (*input->vars)[id+"_s"]=Var(id+"_s", Ek_reduced);
}
