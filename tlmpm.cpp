#include "tlmpm.h"
#include <iostream>
#include <vector>
#include "domain.h"

using namespace std;

TLMPM::TLMPM(MPM *mpm, vector<string> args) : Method(mpm) {
  cout << "In TLMPM::TLMPM()" << endl;
}

TLMPM::~TLMPM(){
  // int nsolids = domain->solids.size();

  // for (int isolid=0; isolid<nsolids; isolid++){
  //   delete neigh_pn[nsolids];
  //   delete neigh_np[nsolids];
  // }
  // delete neigh_pn;
  // delete neigh_np;
}

void TLMPM::setup()
{
  cout << "np = " << domain->solids[0]->np << ", nn = " << domain->solids[0]->grid->nnodes << endl;
  int nsolids, np, nnodes;

  nsolids = domain->solids.size();

  // neigh_pn = new map<int, int>*[nsolids];
  // neigh_np = new map<int, int>*[nsolids];

  // for (int isolid=0; isolid<nsolids; isolid++){
  //   np = domain->solids[isolid]->np;
  //   nnodes = domain->solids[isolid]->grid->nnodes;
  //   neigh_pn[nsolids] = new map<int, int>[np];
  //   neigh_np[nsolids] = new map<int, int>[nnodes];
  // }
}

void TLMPM::compute_grid_weight_functions_and_gradients()
{

}

void TLMPM::particles_to_grid()
{
  /*compute_mass_nodes();
  compute_thermal_energy_nodes();
  compute_velocity_nodes();
  compute_external_forces_nodes();
  compute_internal_forces_nodes();*/
}

void TLMPM::update_grid_state()
{
}

void TLMPM::grid_to_points()
{
}

void TLMPM::advance_particles()
{
}

void TLMPM::velocities_to_grid()
{
}

void TLMPM::compute_rate_deformation_gradient()
{
}

void TLMPM::update_deformation_gradient()
{
}

void TLMPM::update_stress()
{
}
