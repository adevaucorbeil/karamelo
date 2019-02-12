#include "tlmpm.h"
#include "domain.h"
#include "solid.h"
#include "grid.h"
#include "input.h"
#include <iostream>
#include <vector>
#include <Eigen/Eigen>
#include <math.h>

using namespace std;

TLMPM::TLMPM(MPM *mpm, vector<string> args) : Method(mpm) {
  cout << "In TLMPM::TLMPM()" << endl;

  update_wf = 1;
  FLIP = 0.99;
}

TLMPM::~TLMPM()
{
}

void TLMPM::modify(vector<string> args)
{
  FLIP = input->parse(args[0]);
}

void TLMPM::setup()
{
  cout << "np = " << domain->solids[0]->np << ", nn = " << domain->solids[0]->grid->nnodes << endl;
  //compute_grid_weight_functions_and_gradients();
}

void TLMPM::compute_grid_weight_functions_and_gradients()
{
  if (!update_wf) return;

  bigint nsolids, np, nnodes;

  nsolids = domain->solids.size();

  if (nsolids) {
    for (int isolid=0; isolid<nsolids; isolid++){

      np = domain->solids[isolid]->np;
      nnodes = domain->solids[isolid]->grid->nnodes;

      int *numneigh_pn = domain->solids[isolid]->numneigh_pn;
      int *numneigh_np = domain->solids[isolid]->numneigh_np;

      vector<int> *neigh_pn = domain->solids[isolid]->neigh_pn;
      vector<int> *neigh_np = domain->solids[isolid]->neigh_np;

      vector< double > *wf_pn = domain->solids[isolid]->wf_pn;
      vector< double > *wf_np = domain->solids[isolid]->wf_np;

      vector< array<double,3> > *wfd_pn = domain->solids[isolid]->wfd_pn;
      vector< array<double,3> > *wfd_np = domain->solids[isolid]->wfd_np;

      Eigen::Vector3d r;
      double s[3], sd[3];
      Eigen::Vector3d *xp = domain->solids[isolid]->x0;
      Eigen::Vector3d *xn = domain->solids[isolid]->grid->x;
      double inv_cellsize = 1.0 / domain->solids[isolid]->grid->cellsize;
      double wf;
      array<double,3> wfd;

      if (np && nnodes) {
	for (int ip=0; ip<np; ip++) {
	  for (int in=0; in<nnodes; in++) {
	    // Calculate the distance between each pair of particle/node:
	    r = (xp[ip] - xn[in]) * inv_cellsize;

	    s[0] = spline(r[0]);
	    s[1] = spline(r[1]);
	    if (domain->dimension == 3) s[2] = spline(r[2]);

	    if (s[0] != 0 && s[1] != 0 && s[2] != 0) {

	      sd[0] = derivative_spline(r[0], inv_cellsize);
	      sd[1] = derivative_spline(r[1], inv_cellsize);
	      if (domain->dimension == 3) sd[2] = derivative_spline(r[2], inv_cellsize);

	      neigh_pn[ip].push_back(in);
	      neigh_np[in].push_back(ip);
	      numneigh_pn[ip]++;
	      numneigh_np[in]++;

	      if (domain->dimension == 2) wf = s[0]*s[1];
	      if (domain->dimension == 3) wf = s[0]*s[1]*s[2];

	      wf_pn[ip].push_back(wf);
	      wf_np[in].push_back(wf);

	      if (domain->dimension == 2)
		{
		  wfd[0] = sd[0]*s[1];
		  wfd[1] = s[0]*sd[1];
		}
	      else if (domain->dimension == 3)
		{
		  wfd[0] = sd[0]*s[1]*s[2];
		  wfd[1] = s[0]*sd[1]*s[2];
		  wfd[2] = s[0]*s[1]*sd[2];
		}
	      wfd_pn[ip].push_back(wfd);
	      wfd_pn[ip].push_back(wfd);
	    }
	  } 
	}
      }
    }
  }

  update_wf = 0;
}

double TLMPM::spline(double r_)
{
  double r = fabs(r_);
  if (r >= 1.0)
    return 0.0;
  else
    return 1.0 - r;
}

double TLMPM::derivative_spline(double r, double inv_cellsize)
{
  if (r >= 1.0 || r <= -1.0 || r == 0)
    return 0.0;
  else if (r > 0.0)
    return -inv_cellsize;
  else
    return inv_cellsize;
}

void TLMPM::particles_to_grid()
{
  for (int isolid=0; isolid<domain->solids.size(); isolid++){
      domain->solids[isolid]->compute_mass_nodes();
      domain->solids[isolid]->compute_velocity_nodes();
      domain->solids[isolid]->compute_external_forces_nodes();
      domain->solids[isolid]->compute_internal_forces_nodes();
      /*compute_thermal_energy_nodes();*/
    }
}

void TLMPM::update_grid_state()
{
  for (int isolid=0; isolid<domain->solids.size(); isolid++) {
    domain->solids[isolid]->grid->update_grid_velocities();
  }
}

void TLMPM::grid_to_points()
{
  for (int isolid=0; isolid<domain->solids.size(); isolid++) {
    domain->solids[isolid]->compute_particle_velocities();
    domain->solids[isolid]->compute_particle_acceleration();
  }
}

void TLMPM::advance_particles()
{
  for (int isolid=0; isolid<domain->solids.size(); isolid++) {
    domain->solids[isolid]->update_particle_position();
    domain->solids[isolid]->update_particle_velocities(FLIP);
  }
}

void TLMPM::velocities_to_grid()
{
  for (int isolid=0; isolid<domain->solids.size(); isolid++) {
    domain->solids[isolid]->compute_mass_nodes();
    domain->solids[isolid]->compute_velocity_nodes();
  }
}

void TLMPM::compute_rate_deformation_gradient()
{
  for (int isolid=0; isolid<domain->solids.size(); isolid++) {
    domain->solids[isolid]->compute_rate_deformation_gradient();
  }
}

void TLMPM::update_deformation_gradient()
{
  for (int isolid=0; isolid<domain->solids.size(); isolid++) {
    domain->solids[isolid]->update_deformation_gradient();
  }
}

void TLMPM::update_stress()
{
  for (int isolid=0; isolid<domain->solids.size(); isolid++) {
    domain->solids[isolid]->update_stress();
  }
}

