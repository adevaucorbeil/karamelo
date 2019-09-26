#include <iostream>
#include <vector>
#include <Eigen/Eigen>
#include <math.h>
#include <algorithm>
#include "ulcpdi.h"
#include "domain.h"
#include "solid.h"
#include "grid.h"
#include "input.h"
#include "update.h"
#include "var.h"
#include "basis_functions.h"
#include "error.h"

using namespace std;

ULCPDI::ULCPDI(MPM *mpm, vector<string> args) : Method(mpm) {
  cout << "In ULCPDI::ULCPDI()" << endl;

  update_wf = 1;
  method_type = "FLIP";
  FLIP = 0.99;
  is_CPDI = true;

  // Default base function (linear):
  shape_function = "linear";
  basis_function = &BasisFunction::linear;
  derivative_basis_function = &BasisFunction::derivative_linear;
}

ULCPDI::~ULCPDI()
{
}

void ULCPDI::setup(vector<string> args)
{
  int n = 1;
  bool isFLIP = false;
  // Method used: PIC, FLIP or APIC:
  if (args[n].compare("PIC") == 0) {
    method_type = "PIC";
    FLIP = 0;
  } else if (args[n].compare("FLIP") == 0) {
    method_type = "FLIP";
    isFLIP = true;

    if (args.size() < 2) {
      error->all(FLERR, "Illegal modify_method command: not enough arguments.\n");
    }

  } else if (args[n].compare("APIC") == 0) {
    method_type = "APIC";
  } else {
    error->all(FLERR, "Error: method type " + args[n] + " not understood. Expect: PIC, FLIP or APIC\n");
  }

  n++;
  
  if (args.size() > 1 + isFLIP) {
    if (args[n].compare("linear") == 0) {
      shape_function = "linear";
      cout << "Setting up linear basis functions\n";
      basis_function = &BasisFunction::linear;
      derivative_basis_function = &BasisFunction::derivative_linear;
      n++;
    } else if (args[n].compare("cubic-spline") == 0) {
      shape_function = "cubic-spline";
      cout << "Setting up cubic-spline basis functions\n";
      basis_function = &BasisFunction::cubic_spline;
      derivative_basis_function = &BasisFunction::derivative_cubic_spline;
      n++;
    } else if (args[n].compare("Bernstein-quadratic") == 0) {
      shape_function = "Bernstein-quadratic";
      cout << "Setting up Bernstein-quadratic basis functions\n";
      basis_function = &BasisFunction::bernstein_quadratic;
      derivative_basis_function = &BasisFunction::derivative_bernstein_quadratic;
      n++;
    } else {
      error->all(FLERR, "Illegal method_method argument: form function of type " + args[n] + " is unknown.\n");
    }
  }

  if (args.size() > n + isFLIP) {
    error->all(FLERR, "Illegal modify_method command: too many arguments: " + to_string(n + isFLIP) + " expected, " + to_string(args.size()) + " received.\n");
  }

  if (isFLIP) FLIP = input->parsev(args[n]);
  // cout << "shape_function = " << shape_function << endl;
  // cout << "method_type = " << method_type << endl;
  // cout << "FLIP = " << FLIP << endl;
}

void ULCPDI::compute_grid_weight_functions_and_gradients()
{
  if (!update_wf) return;
  if (domain->dimension !=2) {
    error->all(FLERR, "Error: ULCPDI is only 2D....\n");
  }

  bigint nsolids, np, nnodes, nc;

  nsolids = domain->solids.size();

  if (nsolids) {
    for (int isolid=0; isolid<nsolids; isolid++){

      np = domain->solids[isolid]->np;
      nc = domain->solids[isolid]->nc;
      nnodes = domain->solids[isolid]->grid->nnodes;

      int *numneigh_pn = domain->solids[isolid]->numneigh_pn;
      int *numneigh_np = domain->solids[isolid]->numneigh_np;

      vector<int> *neigh_pn = domain->solids[isolid]->neigh_pn;
      vector<int> *neigh_np = domain->solids[isolid]->neigh_np;

      vector< double > *wf_pn = domain->solids[isolid]->wf_pn;
      vector< double > *wf_np = domain->solids[isolid]->wf_np;

      vector< Eigen::Vector3d > *wfd_pn = domain->solids[isolid]->wfd_pn;
      vector< Eigen::Vector3d > *wfd_np = domain->solids[isolid]->wfd_np;

      Eigen::Vector3d *xp = domain->solids[isolid]->x;
      Eigen::Vector3d *xn = domain->solids[isolid]->grid->x0;
      Eigen::Vector3d *rp = domain->solids[isolid]->rp;

      double inv_cellsize = 1.0 / domain->solids[isolid]->grid->cellsize;
      double *vol = domain->solids[isolid]->vol;
      int **ntype = domain->solids[isolid]->grid->ntype;

      double wf;
      double s[3];
      Eigen::Vector3d r, wfd;
      vector<Eigen::Vector3d> xcorner(nc, Eigen::Vector3d::Zero());
      vector<double> wfc(nc, 0);

      for (int in=0; in<nnodes; in++) {
	neigh_np[in].clear();
	numneigh_np[in]=0;
	wf_np[in].clear();
	wfd_np[in].clear();
      }

      if (np && nnodes) {
	for (int ip=0; ip<np; ip++) {

	  neigh_pn[ip].clear();
	  numneigh_pn[ip] = 0;
	  wf_pn[ip].clear();
	  wfd_pn[ip].clear();

	  // Calculate what nodes the corner of Omega_p will interact with:
	  int nx = domain->solids[isolid]->grid->nx;
	  int ny = domain->solids[isolid]->grid->ny;
	  int nz = domain->solids[isolid]->grid->nz;

	  vector<int> n_neigh;
	  int m;


	  // Calculate the coordinates of the particle domain's corners:
	  if (domain->dimension == 1) {
	    xcorner[0] = xp[ip] - rp[ip];
	    xcorner[1] = xp[ip] + rp[ip];
	  }

	  if (domain->dimension == 2) {
	    xcorner[0] = xp[ip] - rp[2*ip] - rp[2*ip+1];
	    xcorner[1] = xp[ip] + rp[2*ip] - rp[2*ip+1];
	    xcorner[2] = xp[ip] + rp[2*ip] + rp[2*ip+1];
	    xcorner[3] = xp[ip] - rp[2*ip] + rp[2*ip+1];
	  }

	  if (domain->dimension == 3) {
	    error->all(FLERR, "Unsupported!\n");
	  }

	  for (int ic=0; ic<nc; ic++) { // Do this for all corners

	    int i0, j0, k0;

	    if (update->method_shape_function.compare("linear")==0) {
	      i0 = (int) ((xcorner[ic][0] - domain->boxlo[0])*inv_cellsize);
	      j0 = (int) ((xcorner[ic][1] - domain->boxlo[1])*inv_cellsize);
	      k0 = (int) ((xcorner[ic][2] - domain->boxlo[2])*inv_cellsize);

	      m = 2;

	    } else if (update->method_shape_function.compare("Bernstein-quadratic")==0){
	      i0 = 2*(int) ((xcorner[ic][0] - domain->boxlo[0])*inv_cellsize);
	      j0 = 2*(int) ((xcorner[ic][1] - domain->boxlo[1])*inv_cellsize);
	      k0 = 2*(int) ((xcorner[ic][2] - domain->boxlo[2])*inv_cellsize);

	      if ((i0 >= 1) && (i0 % 2 != 0)) i0--;
	      if ((j0 >= 1) && (j0 % 2 != 0)) j0--;
	      if (nz>1) if ((k0 >= 1) && (k0 % 2 != 0)) k0--;

	      m = 3;

	    } else if (update->method_shape_function.compare("cubic-spline")==0){
	      i0 = (int) ((xcorner[ic][0] - domain->boxlo[0])*inv_cellsize - 1);
	      j0 = (int) ((xcorner[ic][1] - domain->boxlo[1])*inv_cellsize - 1);
	      k0 = (int) ((xcorner[ic][2] - domain->boxlo[2])*inv_cellsize - 1);

	      m = 4;

	    } else {
	      error->all(FLERR, "Shape function type not supported by TLMPM::compute_grid_weight_functions_and_gradients(): " + update->method_shape_function + ".\n");
	    }

	    // cout << "corner " << ic << " of particle " << ip << ": [" << xcorner[ic][0] << "," << xcorner[ic][1] << "," << xcorner[ic][2] << "], i0=" << i0 << ", j0=" << j0 << ", k0=" << k0 << endl;
	    for(int i=i0; i<i0+m;i++){
	      if (ny>1) {
		for(int j=j0; j<j0+m;j++){
		  if (nz>1){
		    for(int k=k0; k<k0+m;k++){
		      int n = nz*ny*i+nz*j+k;
		      if (n < nnodes)
			n_neigh.push_back(n);
		    }
		  } else {
		    int n = ny*i+j;
		    if (n < nnodes) {
		      // cout << "nodes:" << n << endl;
		      n_neigh.push_back(n);
		    }
		  }
		}
	      } else {
		if (i < nnodes)
		  n_neigh.push_back(i);
	      }
	    }
	  }

	  // Keep only unique values of in in n_neigh:
	  vector<int>::iterator it;
	  sort(n_neigh.begin(), n_neigh.end());           // Sort all values
	  it = unique (n_neigh.begin(), n_neigh.end());   // Remove indentical consecutive values
	  n_neigh.resize( distance(n_neigh.begin(),it) ); // Resize vector.

	  // cout << "[";
	  // for (auto ii: n_neigh)
	  //   cout << ii << ' ';
	  // cout << "]\n";

	  //for (int in=0; in<nnodes; in++) {
	  for (auto in: n_neigh) {
	    wf = 0;

	    for(int ic=0; ic<nc; ic++) {
	      // Calculate the distance between each pair of particle/node:
	      r = (xcorner[ic] - xn[in]) * inv_cellsize;

	      s[0] = basis_function(r[0], ntype[in][0]);
	      s[1] = basis_function(r[1], ntype[in][1]);
	      if (domain->dimension == 3) s[2] = basis_function(r[2], ntype[in][2]);
	      else s[2] = 1;

	      wfc[ic] = s[0]*s[1]*s[2]; // Shape function of the corner node

	      if (wfc[ic] > 1.0e-12) wf += wfc[ic];
	    }

	    wf *= 0.25;

	    if (wf > 1.0e-12) {
	      if (domain->dimension == 2) {
		wfd[0] = (wfc[0] - wfc[2]) * (rp[domain->dimension*ip][1] - rp[domain->dimension*ip+1][1])
		  + (wfc[1] - wfc[3]) * (rp[domain->dimension*ip][1] + rp[domain->dimension*ip+1][1]);
		
		wfd[1] = (wfc[0] - wfc[2]) * (rp[domain->dimension*ip+1][0] - rp[domain->dimension*ip][0])
		  - (wfc[1] - wfc[3]) * (rp[domain->dimension*ip][0] + rp[domain->dimension*ip+1][0]);
		wfd[2] = 0;
	      }

	      double inv_Vp = 1.0/vol[ip];
	      wfd[0] *= inv_Vp;
	      wfd[1] *= inv_Vp;
	      wfd[2] *= inv_Vp;

	      neigh_pn[ip].push_back(in);
	      neigh_np[in].push_back(ip);
	      numneigh_pn[ip]++;
	      numneigh_np[in]++;

	      wf_pn[ip].push_back(wf);
	      wf_np[in].push_back(wf);
	      wfd_pn[ip].push_back(wfd);
	      wfd_np[in].push_back(wfd);
	      // cout << "node: " << in << " [ " << xn[in][0] << "," << xn[in][1] << "," << xn[in][2] << "]" <<
	      // 	" with\twf=" << wf << " and\twfd=["<< wfd[0] << "," << wfd[1] << "," << wfd[2] << "]\n";
	    }
	  }
	}
      }
      if (method_type.compare("APIC") == 0) domain->solids[isolid]->compute_inertia_tensor(shape_function);
    }
  }
}

void ULCPDI::particles_to_grid()
{
  bool grid_reset = false; // Indicate if the grid quantities have to be reset
  for (int isolid=0; isolid<domain->solids.size(); isolid++){

    if (isolid == 0) grid_reset = true;
    else grid_reset = false;

    domain->solids[isolid]->compute_mass_nodes(grid_reset);
    if (method_type.compare("APIC") == 0) domain->solids[isolid]->compute_velocity_nodes_APIC(grid_reset);
    else domain->solids[isolid]->compute_velocity_nodes(grid_reset);
    domain->solids[isolid]->compute_external_forces_nodes(grid_reset);
    domain->solids[isolid]->compute_internal_forces_nodes_UL(grid_reset);
    /*compute_thermal_energy_nodes();*/
    }
}

void ULCPDI::update_grid_state()
{
  domain->grid->update_grid_velocities();
}

void ULCPDI::grid_to_points()
{
  for (int isolid=0; isolid<domain->solids.size(); isolid++) {
    domain->solids[isolid]->compute_particle_velocities();
    domain->solids[isolid]->compute_particle_acceleration();
  }
}

void ULCPDI::advance_particles()
{
  for (int isolid=0; isolid<domain->solids.size(); isolid++) {
    domain->solids[isolid]->update_particle_position();
    domain->solids[isolid]->update_particle_velocities(FLIP);
  }
}

void ULCPDI::velocities_to_grid()
{
  bool grid_reset = false; // Indicate if the grid quantities have to be reset
  for (int isolid=0; isolid<domain->solids.size(); isolid++){

    if (isolid == 0) grid_reset = true;
    else grid_reset = false;

    if (method_type.compare("APIC") != 0) { 
      //domain->solids[isolid]->compute_mass_nodes(grid_reset);
      domain->solids[isolid]->compute_velocity_nodes(grid_reset);
    }
    // domain->solids[isolid]->grid->update_grid_positions();
  }
}

void ULCPDI::compute_rate_deformation_gradient()
{
  for (int isolid=0; isolid<domain->solids.size(); isolid++) {
    if (method_type.compare("APIC") == 0) domain->solids[isolid]->compute_rate_deformation_gradient_UL_APIC();
    else domain->solids[isolid]->compute_rate_deformation_gradient_UL_MUSL();
    //domain->solids[isolid]->compute_deformation_gradient();
  }
}

void ULCPDI::update_deformation_gradient()
{
  for (int isolid=0; isolid<domain->solids.size(); isolid++) {
    domain->solids[isolid]->update_deformation_gradient();
    domain->solids[isolid]->update_particle_domain();
  }
}

void ULCPDI::update_stress()
{
  for (int isolid=0; isolid<domain->solids.size(); isolid++) {
    domain->solids[isolid]->update_stress();
  }
}

void ULCPDI::adjust_dt()
{
  update->update_time();
  if (update->dt_constant) return; // dt is set as a constant, do not update


  double dtCFL = 1.0e22;

  for (int isolid=0; isolid<domain->solids.size(); isolid++) {
    dtCFL = MIN(dtCFL, domain->solids[isolid]->dtCFL);
    if (dtCFL == 0) {
      cout << "Error: dtCFL == 0\n";
      cout << "domain->solids[" << isolid << "]->dtCFL == 0\n";
      error->all(FLERR, "");
    } else if (std::isnan(dtCFL)) {
      cout << "Error: dtCFL = " << dtCFL << "\n";
      cout << "domain->solids[" << isolid << "]->dtCFL == " << domain->solids[isolid]->dtCFL << "\n";
      error->all(FLERR, "");
    }
  }
  update->dt = dtCFL * update->dt_factor;
  (*input->vars)["dt"] = Var("dt", update->dt);
}

void ULCPDI::reset()
{
  int np;

  for (int isolid=0; isolid<domain->solids.size(); isolid++) {
    domain->solids[isolid]->dtCFL = 1.0e22;
    np = domain->solids[isolid]->np;
    for (int ip = 0; ip < np; ip++) domain->solids[isolid]->mb[ip].setZero();
  }
}
