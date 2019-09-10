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
#include <algorithm>
#include "basis_functions.h"

using namespace std;

TLMPM::TLMPM(MPM *mpm, vector<string> args) : Method(mpm) {
  cout << "In TLMPM::TLMPM()" << endl;

  update_wf = 1;
  method_type = "FLIP";
  FLIP = 0.99;

  // Default base function (linear):
  shape_function = "linear";
  basis_function = &BasisFunction::linear;
  derivative_basis_function = &BasisFunction::derivative_linear;
}

TLMPM::~TLMPM()
{
}

void TLMPM::setup(vector<string> args)
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
      cout << "Illegal modify_method command: not enough arguments." << endl;
      exit(1);
    }

  } else if (args[n].compare("APIC") == 0) {
    method_type = "APIC";
  } else {
    cout << "Error: method type " << args[n] << " not understood. Expect: PIC, FLIP or APIC\n";
    exit(1);
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
      cout << "Illegal method_method argument: form function of type " << args[n] << " is unknown." << endl;
      exit(1);
    }
  }

  if (args.size() > n + isFLIP) {
    cout << "Illegal modify_method command: too many arguments: " << n + isFLIP << " expected, " << args.size() << " received." << endl;
      exit(1);    
  }

  if (isFLIP) FLIP = input->parsev(args[n]);
  // cout << "shape_function = " << shape_function << endl;
  // cout << "method_type = " << method_type << endl;
  // cout << "FLIP = " << FLIP << endl;
}

void TLMPM::compute_grid_weight_functions_and_gradients()
{
  if (!update_wf) return;

  cout << "In TLMPM::compute_grid_weight_functions_and_gradients()\n";
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

      int **ntype = domain->solids[isolid]->grid->ntype;

      if (np && nnodes) {
	for (int ip=0; ip<np; ip++) {
	  // Calculate what nodes particle ip will interact with:
	  int nx = domain->solids[isolid]->grid->nx;
	  int ny = domain->solids[isolid]->grid->ny;
	  int nz = domain->solids[isolid]->grid->nz;

	  vector<int> n_neigh;

	  if (update->method_shape_function.compare("linear")==0) {
	    int i0 = (int) ((xp[ip][0] - domain->solids[isolid]->solidlo[0])*inv_cellsize);
	    int j0 = (int) ((xp[ip][1] - domain->solids[isolid]->solidlo[1])*inv_cellsize);
	    int k0 = (int) ((xp[ip][2] - domain->solids[isolid]->solidlo[2])*inv_cellsize);

	    for(int i=i0; i<i0+2;i++){
	      for(int j=j0; j<j0+2;j++){
		if (nz>1){
		  for(int k=k0; k<k0+2;k++){
		    n_neigh.push_back(nz*ny*i+nz*j+k);
		  }
		} else {
		  n_neigh.push_back(ny*i+j);
		}
	      }
	    }
	  } else if (update->method_shape_function.compare("Bernstein-quadratic")==0){
	    int i0 = 2*(int) ((xp[ip][0] - domain->solids[isolid]->solidlo[0])*inv_cellsize);
	    int j0 = 2*(int) ((xp[ip][1] - domain->solids[isolid]->solidlo[1])*inv_cellsize);
	    int k0 = 2*(int) ((xp[ip][2] - domain->solids[isolid]->solidlo[2])*inv_cellsize);

	    if ((i0 >= 1) && (i0 % 2 != 0)) i0--;
	    if ((j0 >= 1) && (j0 % 2 != 0)) j0--;
	    if (nz>1) if ((k0 >= 1) && (k0 % 2 != 0)) k0--;

	    // cout << "(" << i0 << "," << j0 << "," << k0 << ")\t";

	    for(int i=i0; i<i0+3;i++){
	      for(int j=j0; j<j0+3;j++){
		if (nz>1){
		  for(int k=k0; k<k0+3;k++){
		    n_neigh.push_back(nz*ny*i+nz*j+k);
		    // if (nz*ny*i+nz*j+k >= nnodes) {
		    //   cout << "Error: " << nz*ny*i+nz*j+k << " >= nnodes=" << nnodes << endl ;
		    //   exit(1);
		    // }
		  }
		} else {
		  n_neigh.push_back(ny*i+j);
		  // if (ny*i+j >= nnodes) {
		  //   cout << "Error: " << ny*i+j << " >= nnodes=" << nnodes << endl ;
		  //   exit(1);
		  // }
		}
	      }
	    }
	  } else if (update->method_shape_function.compare("cubic-spline")==0){
	    int i0 = (int) ((xp[ip][0] - domain->solids[isolid]->solidlo[0])*inv_cellsize - 1);
	    int j0 = (int) ((xp[ip][1] - domain->solids[isolid]->solidlo[1])*inv_cellsize - 1);
	    int k0 = (int) ((xp[ip][2] - domain->solids[isolid]->solidlo[2])*inv_cellsize - 1);

	    // cout << "(" << i0 << "," << j0 << "," << k0 << ")\t";

	    for(int i=i0; i<i0+4;i++){
	      for(int j=j0; j<j0+4;j++){
		if (nz>1){
		  for(int k=k0; k<k0+4;k++){
		    int n = nz*ny*i+nz*j+k;
		    if (n < nnodes)
		      n_neigh.push_back(n);
		  }
		} else {
		  int n = ny*i+j;
		    if (n < nnodes)
		      n_neigh.push_back(n);
		}
	      }
	    }
	  } else {
	    cout << "Shape function type not supported by TLMPM::compute_grid_weight_functions_and_gradients(): " << update->method_shape_function << endl;
	    exit(1);
	  }

	  // cout << "[";
	  // for (auto ii: n_neigh)
	  //   cout << ii << ' ';
	  // cout << "]\n";

	  //for (int in=0; in<nnodes; in++) {
	  for (auto in: n_neigh) {
	    // Calculate the distance between each pair of particle/node:
	    r = (xp[ip] - xn[in]) * inv_cellsize;

	    s[0] = basis_function(r[0], ntype[in][0]);
	    s[1] = basis_function(r[1], ntype[in][1]);
	    if (domain->dimension == 3) s[2] = basis_function(r[2], ntype[in][2]);
	    else s[2] = 1;

	    if (s[0] != 0 && s[1] != 0 && s[2] != 0) {
	      // // cout << in << "\t";
	      // // Check if this node is in n_neigh:
	      // if (find(n_neigh.begin(), n_neigh.end(), in) == n_neigh.end()) {
	      // 	// in is not in n_neigh
	      //  	cout << "in=" << in << " not found in n_neigh for ip=" << ip << " which is :[";
	      //  	for (auto ii: n_neigh)
	      //  	  cout << ii << ' ';
	      //  	cout << "]\n";
	      // }

	      sd[0] = derivative_basis_function(r[0], ntype[in][0], inv_cellsize);
	      sd[1] = derivative_basis_function(r[1], ntype[in][1], inv_cellsize);
	      if (domain->dimension == 3) sd[2] = derivative_basis_function(r[2], ntype[in][2], inv_cellsize);

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
		  wfd[2] = 0;
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
	  // cout << endl;
	}
      }
      if (method_type.compare("APIC") == 0) domain->solids[isolid]->compute_inertia_tensor(shape_function);
    }
  }

  update_wf = 0;
}

void TLMPM::particles_to_grid()
{
  bool grid_reset = true; // Indicate if the grid quantities have to be reset
  for (int isolid=0; isolid<domain->solids.size(); isolid++){
    domain->solids[isolid]->compute_mass_nodes(grid_reset);
    //domain->solids[isolid]->compute_node_rotation_matrix(grid_reset);
    if (method_type.compare("APIC") == 0) domain->solids[isolid]->compute_velocity_nodes_APIC(grid_reset);
    else domain->solids[isolid]->compute_velocity_nodes(grid_reset);
    domain->solids[isolid]->compute_external_forces_nodes(grid_reset);
    domain->solids[isolid]->compute_internal_forces_nodes_TL();
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
    if (method_type.compare("APIC") != 0) { 
      //domain->solids[isolid]->compute_mass_nodes();
      domain->solids[isolid]->compute_velocity_nodes(true);
    }
    domain->solids[isolid]->grid->update_grid_positions();
  }
}

void TLMPM::compute_rate_deformation_gradient()
{
  for (int isolid=0; isolid<domain->solids.size(); isolid++) {
    if (method_type.compare("APIC") == 0) domain->solids[isolid]->compute_rate_deformation_gradient_TL_APIC();
    else domain->solids[isolid]->compute_rate_deformation_gradient_TL();
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
  update->update_time();
  if (update->dt_constant) return; // dt is set as a constant, do not update


  double dtCFL = 1.0e22;

  for (int isolid=0; isolid<domain->solids.size(); isolid++) {
    dtCFL = MIN(dtCFL, domain->solids[isolid]->dtCFL);
    if (dtCFL == 0) {
      cout << "Error: dtCFL == 0\n";
      cout << "domain->solids[" << isolid << "]->dtCFL == 0\n";
      exit(1);
    } else if (std::isnan(dtCFL)) {
      cout << "Error: dtCFL = " << dtCFL << "\n";
      cout << "domain->solids[" << isolid << "]->dtCFL == " << domain->solids[isolid]->dtCFL << "\n";
      exit(1);
    }
  }
  update->dt = dtCFL * update->dt_factor;
  (*input->vars)["dt"] = Var("dt", update->dt);
}

void TLMPM::reset()
{
  int np;

  for (int isolid=0; isolid<domain->solids.size(); isolid++) {
    domain->solids[isolid]->dtCFL = 1.0e22;
    np = domain->solids[isolid]->np;
    for (int ip = 0; ip < np; ip++) domain->solids[isolid]->b[ip].setZero();
  }
}
