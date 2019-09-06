#include <iostream>
#include <vector>
#include "fix_strain_energy.h"
#include "input.h"
#include "group.h"
#include "domain.h"
#include "input.h"
#include "update.h"
#include "output.h"
#include "math_special.h"
#include <Eigen/Eigen>

using namespace std;
using namespace FixConst;
using namespace MathSpecial;
using namespace Eigen;

FixStrainEnergy::FixStrainEnergy(MPM *mpm, vector<string> args) : Fix(mpm, args)
{
  if (args.size() < 3) {
    cout << "Error: too few arguments for fix_strain_energy: requires at least 3 arguments. " << args.size() << " received" << endl;
    exit(1);
  }

  if (group->pon[igroup].compare("particles") !=0 && group->pon[igroup].compare("all") !=0) {
    cout << "fix_strain_energy needs to be given a group of particles" << group->pon[igroup] << ", " << args[2] << " is a group of "<< group->pon[igroup] << "." << endl;
    exit(1);
  }

  cout << "Creating new fix FixStrainEnergy with ID: " << args[0] << endl;
  id = args[0];
}

FixStrainEnergy::~FixStrainEnergy()
{
}

void FixStrainEnergy::init()
{
}

void FixStrainEnergy::setup()
{
}

void FixStrainEnergy::setmask() {
  mask = 0;
  mask |= FINAL_INTEGRATE;
}


void FixStrainEnergy::final_integrate() {
  if (update->ntimestep != output->next && update->ntimestep != update->nsteps) return;
  // cout << "In FixStrainEnergy::post_particles_to_grid()\n";

  // Go through all the nodes in the group and set b to the right value:
  double ux, uy, uz;
    
  int solid = group->solid[igroup];

  int nmax;
  int *mask;
  double *vol;
  double Es = 0;
  Eigen::Matrix3d *sigma;
  Eigen::Matrix3d *strain;

  if (solid == -1) {
    for (int isolid = 0; isolid < domain->solids.size(); isolid++) {
      vol = domain->solids[isolid]->vol;
      nmax = domain->solids[isolid]->np;
      mask = domain->solids[isolid]->mask;
      sigma = domain->solids[isolid]->sigma;
      strain = domain->solids[isolid]->strain_el;

      for (int in = 0; in < nmax; in++) {
	if (mask[in] & groupbit) {
	  Es += 0.5*vol[in]*(sigma[in](0,0)*strain[in](0,0)
			     + sigma[in](0,1)*strain[in](0,1)
			     + sigma[in](0,2)*strain[in](0,2)
			     + sigma[in](1,0)*strain[in](1,0)
			     + sigma[in](1,1)*strain[in](1,1)
			     + sigma[in](1,2)*strain[in](1,2)
			     + sigma[in](2,0)*strain[in](2,0)
			     + sigma[in](2,1)*strain[in](2,1)
			     + sigma[in](2,2)*strain[in](2,2));
	}
      }
    }


    (*input->vars)[id+"_s"]=Var(id+"_s", Es);
  } else {
    vol = domain->solids[solid]->vol;
    nmax = domain->solids[solid]->np;
    mask = domain->solids[solid]->mask;
    sigma = domain->solids[solid]->sigma;
    strain = domain->solids[solid]->strain_el;

    for (int in = 0; in < nmax; in++) {
      if (mask[in] & groupbit) {
	Es += 0.5*vol[in]*(sigma[in](0,0)*strain[in](0,0)
			   + sigma[in](0,1)*strain[in](0,1)
			   + sigma[in](0,2)*strain[in](0,2)
			   + sigma[in](1,0)*strain[in](1,0)
			   + sigma[in](1,1)*strain[in](1,1)
			   + sigma[in](1,2)*strain[in](1,2)
			   + sigma[in](2,0)*strain[in](2,0)
			   + sigma[in](2,1)*strain[in](2,1)
			   + sigma[in](2,2)*strain[in](2,2));
      }
    }


    (*input->vars)[id+"_s"]=Var(id+"_s", Es);
  }
}
