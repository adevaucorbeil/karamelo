#include "tlmpm.h"
#include "domain.h"
#include "solid.h"
#include "grid.h"
#include "input.h"
#include "update.h"
#include <iostream>
#include <vector>
#include <Eigen/Eigen>
#include <math.h>
#include "var.h"

using namespace std;

TLMPM::TLMPM(MPM *mpm, vector<string> args) : Method(mpm) {
  cout << "In TLMPM::TLMPM()" << endl;

  update_wf = 1;
  FLIP = 0.99;

  // Default base function (linear):
  basis_function = &linear_basis_function;
  derivative_basis_function = &derivative_linear_basis_function;
}

TLMPM::~TLMPM()
{
}

void TLMPM::modify(vector<string> args)
{
  FLIP = input->parsev(args[0]);
  if (args.size() > 1) {
    if (args[1].compare("linear") == 0) {
      cout << "Setting up linear basis functions\n";
      basis_function = &linear_basis_function;
      derivative_basis_function = &derivative_linear_basis_function;
    } else if (args[1].compare("cubic-spline") == 0) {
      cout << "Setting up cubic-spline basis functions\n";
      basis_function = &cubic_spline_basis_function;
      derivative_basis_function = &derivative_cubic_spline_basis_function;
    } else {
      cout << "Illegal run_method argument:" << args[1] << endl;
      exit(1);
    }
  }
}

void TLMPM::setup()
{
  // cout << "np = " << domain->solids[0]->np << ", nn = " << domain->solids[0]->grid->nnodes << endl;
  // compute_grid_weight_functions_and_gradients();
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

      vector< Eigen::Vector3d > *wfd_pn = domain->solids[isolid]->wfd_pn;
      vector< Eigen::Vector3d > *wfd_np = domain->solids[isolid]->wfd_np;

      Eigen::Vector3d r;
      double s[3], sd[3];
      Eigen::Vector3d *xp = domain->solids[isolid]->x0;
      Eigen::Vector3d *xn = domain->solids[isolid]->grid->x0;
      double inv_cellsize = 1.0 / domain->solids[isolid]->grid->cellsize;
      double wf;
      Eigen::Vector3d wfd;

      if (np && nnodes) {
	for (int ip=0; ip<np; ip++) {
	  for (int in=0; in<nnodes; in++) {
	    // Calculate the distance between each pair of particle/node:
	    r = (xp[ip] - xn[in]) * inv_cellsize;

	    s[0] = basis_function(r[0]);
	    s[1] = basis_function(r[1]);
	    if (domain->dimension == 3) s[2] = basis_function(r[2]);

	    if (s[0] != 0 && s[1] != 0 && s[2] != 0) {

	      sd[0] = derivative_basis_function(r[0], inv_cellsize);
	      sd[1] = derivative_basis_function(r[1], inv_cellsize);
	      if (domain->dimension == 3) sd[2] = derivative_basis_function(r[2], inv_cellsize);

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
	      wfd_np[in].push_back(wfd);
	      // cout << "ip=" << ip << ", in=" << in << ", wf=" << wf << ", wfd=[" << wfd[0] << "," << wfd[1] << "," << wfd[2] << "]" << endl;
	    }
	  } 
	}
      }
    }
  }

  update_wf = 0;
}

double linear_basis_function(double r_)
{
  double r = fabs(r_);
  if (r >= 1.0)
    return 0.0;
  else
    return 1.0 - r;
}

double derivative_linear_basis_function(double r, double inv_cellsize)
{
  if (r >= 1.0 || r <= -1.0 || r == 0)
    return 0.0;
  else if (r > 0.0)
    return -inv_cellsize;
  else
    return inv_cellsize;
}

double cubic_spline_basis_function(double r_)
{
  double r = fabs(r_);
  if (r >= 2.0) {
    return 0;
  }

  if (r <= 1.0) {
    return 0.5 * r * r * r - r * r + 2. / 3.;
  } else {
    return -r * r * r / 6. + r * r - 2 * r + 4. / 3.;
  }
}

double derivative_cubic_spline_basis_function(double r_signed, double icellsize)
{
  if (r_signed >= 0.0) {
    
    /*
     * no need to change the sign of r
     */
    
    if (r_signed < 2.0) {
      if (r_signed < 1.0) {
	return icellsize * (1.5 * r_signed * r_signed - 2. * r_signed);
      } else {
	return icellsize * (-0.5 * r_signed * r_signed + 2. * r_signed - 2.);
      }
    } else {
      return 0;
    }
    
  } else {
    double r = fabs(r_signed);
    
    if (r < 2.0) {
      if (r < 1.0) {
	return -icellsize * (1.5 * r * r - 2. * r);
      } else {
	return -icellsize * (-0.5 * r * r + 2. * r - 2.);
      }
    } else {
      return 0;
    }
  }
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
    domain->solids[isolid]->grid->update_grid_positions();
  }
}

void TLMPM::compute_rate_deformation_gradient()
{
  for (int isolid=0; isolid<domain->solids.size(); isolid++) {
    domain->solids[isolid]->compute_rate_deformation_gradient();
    //domain->solids[isolid]->compute_deformation_gradient();
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

void TLMPM::adjust_dt()
{
  if (update->dt_constant) return; // dt is set as a constant, do not update

  update->update_time();

  double dtCFL = 1.0e22;

  for (int isolid=0; isolid<domain->solids.size(); isolid++) {
    dtCFL = MIN(dtCFL, domain->solids[isolid]->dtCFL);
    if (dtCFL == 0) {
      cout << "Error: dtCFL == 0\n";
      cout << "domain->solids[" << isolid << "]->dtCFL == 0\n";
      exit(1);
    } else if (isnan(dtCFL)) {
      cout << "Error: dtCFL = " << dtCFL << "\n";
      cout << "domain->solids[" << isolid << "]->dtCFL == " << domain->solids[isolid]->dtCFL << "\n";
      exit(1);
    }
  }
  update->dt = dtCFL * update->dt_factor;
  (*input->vars)["dt"] = Var("dt", update->dt);
}

void TLMPM::reset_dtCFL()
{
  for (int isolid=0; isolid<domain->solids.size(); isolid++) domain->solids[isolid]->dtCFL = 1.0e22;
}
