#include "fix_max_plastic_strain.h"
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

FixMaxPlasticStrain::FixMaxPlasticStrain(MPM *mpm, vector<string> args)
    : Fix(mpm, args)
{
  if (args.size() < 3)
  {
    cout << "Error: too few arguments for fix_max_plastic_strain: requires at "
            "least "
            "3 arguments. "
         << args.size() << " received" << endl;
    exit(1);
  }

  if (group->pon[igroup].compare("particles") != 0 &&
      group->pon[igroup].compare("all") != 0)
  {
    cout << "fix_max_plastic_strain needs to be given a group of particles"
         << group->pon[igroup] << ", " << args[2] << " is a group of "
         << group->pon[igroup] << "." << endl;
    exit(1);
  }

  cout << "Creating new fix FixMaxPlasticStrain with ID: " << args[0] << endl;
  id = args[0];

  (*input->vars)[id+"_1s"]=Var(id+"_1s", 0);
  (*input->vars)[id+"_2s"]=Var(id+"_2s", 0);
}

FixMaxPlasticStrain::~FixMaxPlasticStrain() {}

void FixMaxPlasticStrain::init() {}

void FixMaxPlasticStrain::setup() {}

void FixMaxPlasticStrain::setmask()
{
  mask = 0;
  mask |= FINAL_INTEGRATE;
}

void FixMaxPlasticStrain::final_integrate()
{
  if (update->ntimestep != output->next && update->ntimestep != update->nsteps)
    return;
  // cout << "In FixMaxPlasticStrain::post_particles_to_grid()\n";

  // Go through all the nodes in the group and set b to the right value:
  double ux, uy, uz;

  int solid = group->solid[igroup];

  int nmax;
  int* mask;
  double* ep;
  double* te;
  double Es(0.), Tmax(0.);

  if (solid == -1)
  {
    for (int isolid = 0; isolid < domain->solids.size(); isolid++)
    {
      ep   = domain->solids[isolid]->eff_plastic_strain;
      te   = domain->solids[isolid]->T;
      nmax = domain->solids[isolid]->np;
      mask = domain->solids[isolid]->mask;

      for (int in = 0; in < nmax; in++)
      {
        if (mask[in] & groupbit)
        {
          Es = max(Es, ep[in]);
          Tmax = max(Tmax, te[in]);
        }
      }
    }

    (*input->vars)[id + "_1s"] = Var(id + "_1s", Es   );
    (*input->vars)[id + "_2s"] = Var(id + "_2s", Tmax );
  }
  else
  {
    ep   = domain->solids[solid]->eff_plastic_strain;
    te   = domain->solids[solid]->T;
    nmax = domain->solids[solid]->np;
    mask = domain->solids[solid]->mask;

    for (int in = 0; in < nmax; in++)
    {
      if (mask[in] & groupbit)
      {
        Es = max(Es, ep[in]);
        Tmax = max(Tmax, te[in]);
      }
    }

    (*input->vars)[id + "_1s"] = Var(id + "_1s", Es   );
    (*input->vars)[id + "_2s"] = Var(id + "_2s", Tmax );
  }
}
