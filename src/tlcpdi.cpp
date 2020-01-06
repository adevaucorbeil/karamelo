#include "tlcpdi.h"
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

TLCPDI::TLCPDI(MPM *mpm, vector<string> args) : Method(mpm) {
  cout << "In TLCPDI::TLCPDI()" << endl;

  update_wf = 1;
  method_type = "FLIP";
  FLIP = 0.99;
  style = 0;    //Default CPDI style is known_styles[style]="R4"; 

  // Default base function (linear):
  shape_function = "linear";
  basis_function = &BasisFunction::linear;
  derivative_basis_function = &BasisFunction::derivative_linear;
}

TLCPDI::~TLCPDI()
{
}

void TLCPDI::setup(vector<string> args)
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

  if (args.size() > n + isFLIP + 1) {
    cout << "Illegal modify_method command: too many arguments: " << n + isFLIP + 1 << " expected, " << args.size() << " received." << endl;
      exit(1);    
  }

  if (isFLIP) FLIP = input->parsev(args[n]);
  // cout << "shape_function = " << shape_function << endl;
  // cout << "method_type = " << method_type << endl;
  // cout << "FLIP = " << FLIP << endl;
  
  n++;

  if (n<args.size()) {
    bool found_style = false;
    for (int i=0; i<sizeof(known_styles)/sizeof(string); i++) {
      if (known_styles[i].compare(args[n])==0) {
        style = i;
  found_style = true;
        break;
      }
    }
    if (!found_style) {
      cout << "CPDI style \033[1;31m" << args[n] << "\033[0m unknown. Available options are:";
      for (int i=0; i<sizeof(known_styles)/sizeof(string); i++) {
  if (i) cout << ",";
  cout << " \033[1;32m" << known_styles[i] << "\033[0m";
      }
      cout << ".\n";
      exit(1);
    }
  }

  cout << "Using CPDI-" << known_styles[style] << endl;
}

void TLCPDI::compute_grid_weight_functions_and_gradients()
{
  if (!update_wf) return;

  if (domain->dimension !=2) {
    cout << "Error: ULCPDI is only 2D....\n";
    exit(1);
  }
  cout << "In TLCPDI::compute_grid_weight_functions_and_gradients()\n";
  bigint nsolids, np_local, nnodes, nc;

  nsolids = domain->solids.size();

  if (nsolids) {
    for (int isolid=0; isolid<nsolids; isolid++){

      np_local = domain->solids[isolid]->np_local;
      nc = domain->solids[isolid]->nc;
      nnodes = domain->solids[isolid]->grid->nnodes;

      vector<int> *numneigh_pn = &domain->solids[isolid]->numneigh_pn;
      vector<int> *numneigh_np = &domain->solids[isolid]->numneigh_np;

      vector<vector<int>> *neigh_pn = &domain->solids[isolid]->neigh_pn;
      vector<vector<int>> *neigh_np = &domain->solids[isolid]->neigh_np;

      vector<vector< double >> *wf_pn = &domain->solids[isolid]->wf_pn;
      vector<vector< double >> *wf_pn_corners = &domain->solids[isolid]->wf_pn_corners;
      vector<vector< double >> *wf_np = &domain->solids[isolid]->wf_np;

      vector<vector< Eigen::Vector3d >> *wfd_pn = &domain->solids[isolid]->wfd_pn;
      vector<vector< Eigen::Vector3d >> *wfd_np = &domain->solids[isolid]->wfd_np;

      vector<Eigen::Vector3d> *xp  = &domain->solids[isolid]->x0;
      vector<Eigen::Vector3d> *xpc = &domain->solids[isolid]->xpc;
      vector<Eigen::Vector3d> *xn  = &domain->solids[isolid]->grid->x0;
      vector<Eigen::Vector3d> *rp  = &domain->solids[isolid]->rp;

      double inv_cellsize          = 1.0 / domain->solids[isolid]->grid->cellsize;
      vector<double> *vol          = &domain->solids[isolid]->vol;
      vector<array<int, 3>> *ntype = &domain->solids[isolid]->grid->ntype;

      double wf;
      double s[3];
      Eigen::Vector3d r, wfd;
      vector<Eigen::Vector3d> xcorner(nc, Eigen::Vector3d::Zero());
      vector<double> wfc(nc, 0);

      bool linear, cubic, bernstein;
      linear = cubic = bernstein = false;

      if (update->method_shape_function.compare("linear")==0) linear = true;
      if (update->method_shape_function.compare("cubic-spline")==0) cubic = true;
      if (update->method_shape_function.compare("Bernstein-quadratic")==0) bernstein = true;

      double a, b, inv_Vp, alpha_over_Vp, sixVp;

      for (int in=0; in<nnodes; in++) {
	(*neigh_np)[in].clear();
	(*numneigh_np)[in] = 0;
	(*wf_np)[in].clear();
	(*wfd_np)[in].clear();
      }

      if (np_local && nnodes) {
	for (int ip=0; ip<np_local; ip++) {

	  (*neigh_pn)[ip].clear();
	  (*numneigh_pn)[ip] = 0;
	  (*wf_pn)[ip].clear();
	  (*wfd_pn)[ip].clear();

	  // Calculate what nodes the corner of Omega_p will interact with:
	  int nx = domain->solids[isolid]->grid->nx;
	  int ny = domain->solids[isolid]->grid->ny;
	  int nz = domain->solids[isolid]->grid->nz;

	  vector<int> n_neigh;
	  int m;

	  if (style==0) { //CPDI-R4
	    // Calculate the coordinates of the particle domain's corners:
	    if (domain->dimension == 1) {
	      xcorner[0] = (*xp)[ip] - (*rp)[ip];
	      xcorner[1] = (*xp)[ip] + (*rp)[ip];
	    }

	    if (domain->dimension == 2) {
	      xcorner[0] = (*xp)[ip] - (*rp)[2*ip] - (*rp)[2*ip+1];
	      xcorner[1] = (*xp)[ip] + (*rp)[2*ip] - (*rp)[2*ip+1];
	      xcorner[2] = (*xp)[ip] + (*rp)[2*ip] + (*rp)[2*ip+1];
	      xcorner[3] = (*xp)[ip] - (*rp)[2*ip] + (*rp)[2*ip+1];
	    }

	    if (domain->dimension == 3) {
	      cout << "Unsupported!\n";
	      exit(1);
	    }
	  }

	  for (int ic=0; ic<nc; ic++) { // Do this for all corners

	    if (style==1) { // CPDI-Q4
	      xcorner[ic][0] = (*xpc)[nc*ip+ic][0];
	      xcorner[ic][1] = (*xpc)[nc*ip+ic][1];
	      xcorner[ic][2] = (*xpc)[nc*ip+ic][2];
	    }

	    int i0, j0, k0;

	    if (linear) {
	      i0 = (int) ((xcorner[ic][0] - domain->boxlo[0])*inv_cellsize);
	      j0 = (int) ((xcorner[ic][1] - domain->boxlo[1])*inv_cellsize);
	      k0 = (int) ((xcorner[ic][2] - domain->boxlo[2])*inv_cellsize);

	      m = 2;

	    } else if (bernstein){
	      i0 = 2*(int) ((xcorner[ic][0] - domain->boxlo[0])*inv_cellsize);
	      j0 = 2*(int) ((xcorner[ic][1] - domain->boxlo[1])*inv_cellsize);
	      k0 = 2*(int) ((xcorner[ic][2] - domain->boxlo[2])*inv_cellsize);

	      if ((i0 >= 1) && (i0 % 2 != 0)) i0--;
	      if ((j0 >= 1) && (j0 % 2 != 0)) j0--;
	      if (nz>1) if ((k0 >= 1) && (k0 % 2 != 0)) k0--;

	      m = 3;

	    } else if (cubic){
	      i0 = (int) ((xcorner[ic][0] - domain->boxlo[0])*inv_cellsize - 1);
	      j0 = (int) ((xcorner[ic][1] - domain->boxlo[1])*inv_cellsize - 1);
	      k0 = (int) ((xcorner[ic][2] - domain->boxlo[2])*inv_cellsize - 1);

	      m = 4;

	    } else {
	      cout << "Shape function type not supported by ULCPDI::compute_grid_weight_functions_and_gradients(): " << update->method_shape_function << endl;
	      exit(1);
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

	  inv_Vp = 1.0 / (*vol)[ip];

	  if (style==1) {
	    a = (xcorner[3][0]-xcorner[0][0])*(xcorner[1][1]-xcorner[2][1])
	      - (xcorner[1][0]-xcorner[2][0])*(xcorner[3][1]-xcorner[0][1]);

	    b = (xcorner[2][0]-xcorner[3][0])*(xcorner[0][1]-xcorner[1][1])
	      - (xcorner[0][0]-xcorner[1][0])*(xcorner[2][1]-xcorner[3][1]);

	    alpha_over_Vp = 0.0417 * inv_Vp;
	    sixVp = 6 * (*vol)[ip];
	  }

	  //for (int in=0; in<nnodes; in++) {
	  for (auto in: n_neigh) {
	    wf = 0;

	    for(int ic=0; ic<nc; ic++) {
	      // Calculate the distance between each pair of particle/node:
	      r = (xcorner[ic] - (*xn)[in]) * inv_cellsize;

	      s[0] = basis_function(r[0], (*ntype)[in][0]);
	      s[1] = basis_function(r[1], (*ntype)[in][1]);
	      if (domain->dimension == 3) s[2] = basis_function(r[2], (*ntype)[in][2]);
	      else s[2] = 1;

	      wfc[ic] = s[0]*s[1]*s[2]; // Shape function of the corner node

	      if (style==0 && wfc[ic] > 1.0e-12) wf += wfc[ic];
	    }

	    if (style==0) wf *= 0.25;
	    if (style==1) {
	      wf = alpha_over_Vp*((sixVp - a - b)*wfc[0]
				  + (sixVp - a + b)*wfc[1]
				  + (sixVp + a + b)*wfc[2]
				  + (sixVp + a - b)*wfc[3]);
	    }

	    
	    if (wf > 1.0e-12) {
	      if (style==0) {
		wfd[0] = (wfc[0] - wfc[2]) * ((*rp)[domain->dimension*ip][1] - (*rp)[domain->dimension*ip+1][1])
		  + (wfc[1] - wfc[3]) * ((*rp)[domain->dimension*ip][1] + (*rp)[domain->dimension*ip+1][1]);
		
		wfd[1] = (wfc[0] - wfc[2]) * ((*rp)[domain->dimension*ip+1][0] - (*rp)[domain->dimension*ip][0])
		  - (wfc[1] - wfc[3]) * ((*rp)[domain->dimension*ip][0] + (*rp)[domain->dimension*ip+1][0]);

		wfd[2] = 0;

		wfd *= inv_Vp;
	      }

	      if (style==1) {
		wfd[0] = wfc[0] * (xcorner[1][1]-xcorner[3][1])
		  + wfc[1] * (xcorner[2][1]-xcorner[0][1])
		  + wfc[2] * (xcorner[3][1]-xcorner[1][1])
		  + wfc[3] * (xcorner[0][1]-xcorner[2][1]);
		
		wfd[1] = wfc[0] * (xcorner[3][0]-xcorner[1][0])
		  + wfc[1] * (xcorner[0][0]-xcorner[2][0])
		  + wfc[2] * (xcorner[1][0]-xcorner[3][0])
		  + wfc[3] * (xcorner[2][0]-xcorner[0][0]);

		wfd[2] = 0;

		wfd *= 0.5*inv_Vp;
		for(int ic=0; ic<nc; ic++) (*wf_pn_corners)[nc*ip+ic].push_back(wfc[ic]);
	      }

	      (*neigh_pn)[ip].push_back(in);
	      (*neigh_np)[in].push_back(ip);
	      (*numneigh_pn)[ip]++;
	      (*numneigh_np)[in]++;

	      (*wf_pn)[ip].push_back(wf);
	      (*wf_np)[in].push_back(wf);
	      (*wfd_pn)[ip].push_back(wfd);
	      (*wfd_np)[in].push_back(wfd);
	      // cout << "node: " << in << " [ " << (*xn)[in][0] << "," << (*xn)[in][1] << "," << (*xn)[in][2] << "]" <<
	      // 	" with\twf=" << wf << " and\twfd=["<< wfd[0] << "," << wfd[1] << "," << wfd[2] << "]\n";
	    }
	  }
	}
      }
      if (method_type.compare("APIC") == 0) domain->solids[isolid]->compute_inertia_tensor(shape_function);
    }
  }
  update_wf = 0;
}

void TLCPDI::particles_to_grid()
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

void TLCPDI::update_grid_state()
{
  for (int isolid=0; isolid<domain->solids.size(); isolid++) {
    domain->solids[isolid]->grid->update_grid_velocities();
  }
}

void TLCPDI::grid_to_points()
{
  for (int isolid=0; isolid<domain->solids.size(); isolid++) {
    domain->solids[isolid]->compute_particle_velocities_and_positions();
    domain->solids[isolid]->compute_particle_acceleration();
  }
}

void TLCPDI::advance_particles()
{
  for (int isolid=0; isolid<domain->solids.size(); isolid++) {
    domain->solids[isolid]->update_particle_velocities(FLIP);
  }
}

void TLCPDI::velocities_to_grid()
{
  if (method_type.compare("APIC") != 0) {
    for (int isolid=0; isolid<domain->solids.size(); isolid++) {
      //domain->solids[isolid]->compute_mass_nodes();
      domain->solids[isolid]->compute_velocity_nodes(true);
    }
  }
}

void TLCPDI::update_grid_positions()
{
  for (int isolid=0; isolid<domain->solids.size(); isolid++) {
    domain->solids[isolid]->grid->update_grid_positions();
  }
}

void TLCPDI::compute_rate_deformation_gradient()
{
  for (int isolid=0; isolid<domain->solids.size(); isolid++) {
    if (method_type.compare("APIC") == 0) domain->solids[isolid]->compute_rate_deformation_gradient_TL_APIC();
    else domain->solids[isolid]->compute_rate_deformation_gradient_TL();
    //domain->solids[isolid]->compute_deformation_gradient();
  }
}

void TLCPDI::update_deformation_gradient()
{
  for (int isolid=0; isolid<domain->solids.size(); isolid++) {
    domain->solids[isolid]->update_deformation_gradient();
  }
}

void TLCPDI::update_stress()
{
  for (int isolid=0; isolid<domain->solids.size(); isolid++) {
    domain->solids[isolid]->update_stress();
  }
}

void TLCPDI::adjust_dt()
{
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

void TLCPDI::reset()
{
  int np_local;

  for (int isolid=0; isolid<domain->solids.size(); isolid++) {
    domain->solids[isolid]->dtCFL = 1.0e22;
    np_local = domain->solids[isolid]->np_local;
    for (int ip = 0; ip < np_local; ip++) domain->solids[isolid]->mbp[ip].setZero();
  }
}
