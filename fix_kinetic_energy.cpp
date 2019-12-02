#include "fix_kinetic_energy.h"
#include "domain.h"
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

FixKineticEnergy::FixKineticEnergy(MPM *mpm, vector<string> args)
    : Fix(mpm, args)
{
  if (args.size() < 3)
  {
    cout << "Error: too few arguments for fix_kinetic_energy: requires at "
            "least 3 arguments. "
         << args.size() << " received" << endl;
    exit(1);
  }

  if (group->pon[igroup].compare("particles") != 0 &&
      group->pon[igroup].compare("all") != 0)
  {
    cout << "fix_kinetic_energy needs to be given a group of particles"
         << group->pon[igroup] << ", " << args[2] << " is a group of "
         << group->pon[igroup] << "." << endl;
    exit(1);
  }

  cout << "Creating new fix FixKineticEnergy with ID: " << args[0] << endl;
  id = args[0];
}

FixKineticEnergy::~FixKineticEnergy() {}

void FixKineticEnergy::init() {}

void FixKineticEnergy::setup() {}

void FixKineticEnergy::setmask()
{
  mask = 0;
  mask |= FINAL_INTEGRATE;
}

void FixKineticEnergy::final_integrate()
{
  if (update->ntimestep != output->next && update->ntimestep != update->nsteps)
    return;
  // cout << "In FixKineticEnergy::post_particles_to_grid()\n";

  // Go through all the nodes in the group and set b to the right value:
  double ux, uy, uz;

  int solid = group->solid[igroup];

  int nmax;
  int *mask;
  double *mass;
  double Ek = 0;
  Eigen::Vector3d *v;

  if (solid == -1)
  {
    for (int isolid = 0; isolid < domain->solids.size(); isolid++)
    {
      mass = domain->solids[isolid]->mass;
      nmax = domain->solids[isolid]->np;
      mask = domain->solids[isolid]->mask;
      v    = domain->solids[isolid]->v;

      for (int in = 0; in < nmax; in++)
      {
        if (mask[in] & groupbit)
        {
          Ek += 0.5 * mass[in] * square(v[in].norm());
          // cout << "in=" << in << "\tmass=" << mass[in] << "\tv[in]=[" <<
          // v[in](0) << "," << v[in](1) << "," << v[in](2) << "]\tdEk=" <<
          // 0.5*mass[in]*v[in].norm() << endl;
        }
      }
    }

    (*input->vars)[id + "_s"] = Var(id + "_s", Ek);
  }
  else
  {
    mass = domain->solids[solid]->mass;
    nmax = domain->solids[solid]->np;
    mask = domain->solids[solid]->mask;
    v    = domain->solids[solid]->v;

    for (int in = 0; in < nmax; in++)
    {
      if (mask[in] & groupbit)
      {
        Ek += 0.5 * mass[in] * square(v[in].norm());
      }
    }

    (*input->vars)[id + "_s"] = Var(id + "_s", Ek);
  }
}
